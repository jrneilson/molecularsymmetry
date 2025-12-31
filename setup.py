#!/usr/bin/env python3
"""
Setup script for molecular symmetry analysis package.
"""

from setuptools import setup, find_packages
import os

# Read the README file for the long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements if they exist
requirements = []
if os.path.exists("requirements.txt"):
    with open("requirements.txt", "r", encoding="utf-8") as fh:
        requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="molecular-symmetry",
    version="1.1.0",
    author="J.R. Neilson",
    author_email="james.neilson@colostate.edu",  # Update with your email
    description="A comprehensive Python package for character table analysis and molecular orbital symmetry reduction using Schoenflies notation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jrneilson/molecularsymmetry",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
        ],
        "jupyter": [
            "jupyter",
            "matplotlib",
        ],
    },
    keywords="chemistry molecular-orbital symmetry character-table point-group schoenflies",
    project_urls={
        "Bug Reports": "https://github.com/jrneilson/molecularsymmetry/issues",
        "Source": "https://github.com/jrneilson/molecularsymmetry",
        "Documentation": "https://github.com/jrneilson/molecularsymmetry#readme",
    },
    entry_points={
        "console_scripts": [
            "molecular-symmetry-demo=molecular_symmetry.examples:main",
        ],
    },
)