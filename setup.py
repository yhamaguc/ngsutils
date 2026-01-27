import sys
import os
import re
import ftplib

from setuptools import setup, find_packages


setup(
    packages=find_packages(),
    package_data={
        'ngsutils': ['data/*.pkl'],
    },
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'extract_upper_lower=ngsutils.extract_upper_lower:main',
            'gtf2sqlite=ngsutils.gtf2sqlite:main',
            'gtf2bed4igv=ngsutils.gtf2bed4igv:main',
            'sqlite2gtf=ngsutils.sqlite2gtf:main',
            'gtf2tsv=ngsutils.gtf2tsv:main',
            'gene2tx=ngsutils.gene2tx:main',
            'id2name=ngsutils.id2name:main',
            'name2id=ngsutils.name2id:main',
            'parse_gtf=ngsutils.parse_gtf:main',
            'remove_chr_from_reference.py=ngsutils.remove_chr_from_reference:main',
            'maf2bed=ngsutils.maf2bed:main',
            'format_transcripts_fasta=ngsutils.format_transcripts_fasta:main',
            'format_refseq_gff.py=format_refseq_gff:main',
            'append_biomart_to_annotation_sqlite=ngsutils.append_biomart_to_annotation_sqlite:main',
            'rm_tagfrombam=ngsutils.rm_tagfrombam:main',
            'extract_splice_sites=ngsutils.extract_splice_sites:main',
        ]
    },
)
