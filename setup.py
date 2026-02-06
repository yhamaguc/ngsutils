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
)
