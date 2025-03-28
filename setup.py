#!/usr/bin/env python
# encoding: utf-8

from setuptools import find_packages, setup
import os

setup(
    name="mcp_us_airline",
    version="0.1",
    description="minimum-cost-percolation model for US airline networks",
    author="Minsuk Kim",
    package_dir={"": "libs/mcp_us_airline"},
    packages=find_packages(where="libs/mcp_us_airline"),
)