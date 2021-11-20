#!/usr/bin/env python
import sys

from setuptools import setup
from setuptools_rust import RustExtension

setup(
    rust_extensions=[RustExtension("rust_ext.rust_ext")],
)
