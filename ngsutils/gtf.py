# Ported from gtfparse (https://github.com/openvax/gtfparse)
# Licensed under the Apache License, Version 2.0
#
# Original authors: Alex Rubinsteyn et al.
# Ported into ngsutils to remove the gtfparse dependency,
# which did not support numpy > 2.

import logging
from collections import OrderedDict
from os.path import exists
from sys import intern

import polars
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

__version__ = "2.5.0"


# ---------------------------------------------------------------------------
# ParsingError
# ---------------------------------------------------------------------------

class ParsingError(Exception):
    pass


# ---------------------------------------------------------------------------
# attribute_parsing
# ---------------------------------------------------------------------------

def expand_attribute_strings(
        attribute_strings,
        quote_char="'",
        missing_value="",
        usecols=None):
    """
    The last column of GTF has a variable number of key value pairs
    of the format: "key1 value1; key2 value2;"
    Parse these into a dictionary mapping each key onto a list of values,
    where the value is None for any row where the key was missing.
    """
    n = len(attribute_strings)

    extra_columns = {}
    column_order = []

    column_interned_strings = {}

    for (i, kv_strings) in enumerate(attribute_strings):
        if type(kv_strings) is str:
            kv_strings = kv_strings.split(";")
        for kv in kv_strings:
            parts = kv.strip().split(" ", 2)[:2]

            if len(parts) != 2:
                continue

            column_name, value = parts

            try:
                column_name = column_interned_strings[column_name]
            except KeyError:
                column_name = intern(str(column_name))
                column_interned_strings[column_name] = column_name

            if usecols is not None and column_name not in usecols:
                continue

            if value[0] == quote_char:
                value = value.replace(quote_char, "")

            try:
                column = extra_columns[column_name]
                old_value = column[i]
                if old_value is missing_value:
                    column[i] = value
                else:
                    column[i] = "%s,%s" % (old_value, value)
            except KeyError:
                column = [missing_value] * n
                column[i] = value
                extra_columns[column_name] = column
                column_order.append(column_name)

    logging.info("Extracted GTF attributes: %s" % column_order)
    return OrderedDict(
        (column_name, extra_columns[column_name])
        for column_name in column_order)


# ---------------------------------------------------------------------------
# create_missing_features
# ---------------------------------------------------------------------------

def create_missing_features(
        dataframe,
        unique_keys={},
        extra_columns={},
        missing_value=None):
    """
    Helper function used to construct a missing feature such as 'transcript'
    or 'gene'. Some GTF files only have 'exon' and 'CDS' entries, but have
    transcript_id and gene_id annotations which allow us to construct those
    missing features.
    """
    if hasattr(dataframe, "to_pandas"):
        dataframe = dataframe.to_pandas()

    extra_dataframes = []

    existing_features = set(dataframe["feature"])
    existing_columns = set(dataframe.columns)

    for (feature_name, groupby_key) in unique_keys.items():

        if feature_name in existing_features:
            logging.info(
                "Feature '%s' already exists in GTF data" % feature_name)
            continue
        logging.info("Creating rows for missing feature '%s'" % feature_name)

        missing = pd.Series([
            x is None or x == ""
            for x in dataframe[groupby_key]])
        not_missing = ~missing
        row_groups = dataframe[not_missing].groupby(groupby_key)

        feature_values = OrderedDict([
            (column_name, [missing_value] * row_groups.ngroups)
            for column_name in dataframe.keys()
        ])

        feature_columns = list(extra_columns.get(feature_name, []))

        for i, (feature_id, group) in enumerate(row_groups):
            feature_values["feature"][i] = feature_name
            feature_values[groupby_key][i] = feature_id
            feature_values["source"][i] = "gtfparse"
            feature_values["start"][i] = group["start"].min()
            feature_values["end"][i] = group["end"].max()
            feature_values["seqname"][i] = group["seqname"].iat[0]
            feature_values["strand"][i] = group["strand"].iat[0]

            for column_name in feature_columns:
                if column_name not in existing_columns:
                    raise ValueError(
                        "Column '%s' does not exist in GTF, columns = %s" % (
                            column_name, existing_columns))

                unique_values = group[column_name].dropna().unique()
                if len(unique_values) == 1:
                    feature_values[column_name][i] = unique_values[0]
        extra_dataframes.append(pd.DataFrame(feature_values))
    return pd.concat([dataframe] + extra_dataframes, ignore_index=True)


# ---------------------------------------------------------------------------
# read_gtf (core)
# ---------------------------------------------------------------------------

REQUIRED_COLUMNS = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]

DEFAULT_COLUMN_DTYPES = {
    "seqname": polars.Categorical,
    "source": polars.Categorical,
    "start": polars.Int64,
    "end": polars.Int64,
    "score": polars.Float32,
    "feature": polars.Categorical,
    "strand": polars.Categorical,
    "frame": polars.UInt32,
}


def parse_with_polars_lazy(
        filepath_or_buffer,
        split_attributes=True,
        features=None,
        fix_quotes_columns=["attribute"]):
    polars.enable_string_cache()
    kwargs = dict(
        has_header=False,
        separator="\t",
        comment_prefix="#",
        null_values=".",
        schema_overrides=DEFAULT_COLUMN_DTYPES)
    try:
        df = polars.read_csv(
                filepath_or_buffer,
                new_columns=REQUIRED_COLUMNS,
                **kwargs).lazy()
    except polars.exceptions.ShapeError:
        raise ParsingError("Wrong number of columns")

    df = df.with_columns([
        polars.col("frame").fill_null(0),
        polars.col("attribute").str.replace_all('"', "'")
    ])

    for fix_quotes_column in fix_quotes_columns:
        df = df.with_columns([
            polars.col(fix_quotes_column).str.replace(';\"', '\"').str.replace(";-", "-")
        ])

    if features is not None:
        features = sorted(set(features))
        df = df.filter(polars.col("feature").is_in(features))

    if split_attributes:
        df = df.with_columns([
            polars.col("attribute").str.split(";").alias("attribute_split")
        ])
    return df


def parse_gtf(
        filepath_or_buffer,
        split_attributes=True,
        features=None,
        fix_quotes_columns=["attribute"]):
    df_lazy = parse_with_polars_lazy(
        filepath_or_buffer=filepath_or_buffer,
        split_attributes=split_attributes,
        features=features,
        fix_quotes_columns=fix_quotes_columns)
    return df_lazy.collect()


def parse_gtf_pandas(*args, **kwargs):
    return parse_gtf(*args, **kwargs).to_pandas()


def parse_gtf_and_expand_attributes(
        filepath_or_buffer,
        restrict_attribute_columns=None,
        features=None):
    """
    Parse lines into column->values dictionary and then expand
    the 'attribute' column into multiple columns.
    """
    df = parse_gtf(
        filepath_or_buffer=filepath_or_buffer,
        features=features,
        split_attributes=True)
    if type(restrict_attribute_columns) is str:
        restrict_attribute_columns = {restrict_attribute_columns}
    elif restrict_attribute_columns:
        restrict_attribute_columns = set(restrict_attribute_columns)
    df.drop_in_place("attribute")
    attribute_pairs = df.drop_in_place("attribute_split")
    return df.with_columns([
        polars.Series(k, vs)
        for (k, vs) in
        expand_attribute_strings(attribute_pairs).items()
        if restrict_attribute_columns is None or k in restrict_attribute_columns
    ])


def read_gtf(
        filepath_or_buffer,
        expand_attribute_column=True,
        infer_biotype_column=False,
        column_converters={},
        column_cast_types={},
        usecols=None,
        features=None,
        result_type='polars'):
    """
    Parse a GTF into a Polars DataFrame (default), Pandas DataFrame, or dict.

    Parameters
    ----------
    filepath_or_buffer : str or buffer object
        Path to GTF file (may be gzip compressed) or buffer object.

    expand_attribute_column : bool
        Replace the 'attribute' column with one column per distinct key.

    infer_biotype_column : bool
        If the 'source' column contains biotype values (e.g. 'protein_coding'),
        copy it to gene_biotype / transcript_biotype columns.

    column_converters : dict, optional
        Mapping from column names to conversion functions.

    column_cast_types : dict, optional
        Mapping from column names to Polars dtypes.

    usecols : list of str or None
        Restrict which columns are returned.

    features : set of str or None
        Drop rows whose feature type is not in this set.

    result_type : 'polars', 'pandas', or 'dict'
        Output format. Default is 'polars'.
    """
    if type(filepath_or_buffer) is str and not exists(filepath_or_buffer):
        raise ValueError("GTF file does not exist: %s" % filepath_or_buffer)

    if expand_attribute_column:
        result_df = parse_gtf_and_expand_attributes(
            filepath_or_buffer,
            restrict_attribute_columns=usecols,
            features=features)
    else:
        result_df = parse_gtf(filepath_or_buffer, features=features)

    result_df = result_df.to_pandas()

    if column_converters or column_cast_types:
        def wrap_to_always_accept_none(f):
            def wrapped_fn(x):
                if x is None or x == "":
                    return None
                else:
                    return f(x)
            return wrapped_fn

        column_names = set(column_converters.keys()).union(column_cast_types.keys())
        for column_name in column_names:
            if column_name in column_converters:
                column_fn = wrap_to_always_accept_none(
                    column_converters[column_name])
                result_df[column_name] = result_df[column_name].apply(column_fn)

            if column_name in column_cast_types:
                column_type = column_cast_types[column_name]
                result_df[column_name] = result_df[column_name].astype(column_type)

    if infer_biotype_column:
        unique_source_values = set(result_df["source"])
        if "protein_coding" in unique_source_values:
            column_names = set(result_df.columns)
            if "gene_biotype" not in column_names:
                logging.info("Using column 'source' to replace missing 'gene_biotype'")
                result_df['gene_biotype'] = result_df['source']
            if "transcript_biotype" not in column_names:
                logging.info("Using column 'source' to replace missing 'transcript_biotype'")
                result_df['transcript_biotype'] = result_df['source']

    if usecols is not None:
        column_names = set(result_df.columns)
        valid_columns = [c for c in usecols if c in column_names]
        result_df = result_df[valid_columns]

    if result_type == "pandas":
        return result_df
    elif result_type == "polars":
        return polars.from_pandas(result_df)
    elif result_type == "dict":
        return result_df.to_dict()
    return result_df
