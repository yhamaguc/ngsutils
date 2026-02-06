"""
ngsutils.utils
~~~~~~~~~~~~~~

This module provides utility functions

"""

import os

import yaml


def root_path():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def from_root(relpath):
    return os.path.join(root_path(), relpath)


def load_conf(relpath):
    with open(from_root(relpath)) as f:
        dict_yaml = yaml.load(f)
    return dict_yaml
