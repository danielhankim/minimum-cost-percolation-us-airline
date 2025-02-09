#!/usr/bin/env python
# encoding: utf-8

from setuptools import find_packages, setup

setup(
    name="mcp_us_airline",
    version="0.1",
    description="minimum-cost-percolation model for US airline networks",
    author="Minsuk Kim",
    packages=find_packages(exclude=("tests",)),
    install_requires=[
        "numpy",
        "pandas",
        "scikit-learn",
        "geopy",
        "geopandas",
        "scipy",
        "tqdm",
        "rasterio",
        "timezonefinder",
        "pytz",
        "shapely"
    ]
)
