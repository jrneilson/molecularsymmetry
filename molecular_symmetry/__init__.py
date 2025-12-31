"""
Molecular Symmetry Analysis Package

A comprehensive Python package for character table analysis and molecular orbital 
symmetry reduction using Schoenflies notation.

Author: J.R. Neilson
Version: 1.1.0
"""

from .molecular_symmetry import (
    PointGroup,
    PointGroupFactory,
    get_point_group,
    # Point group classes
    C1, Ci, Cs, C2, C3, C2v, C3v, C4v, C5v, C6v,
    D2, D3, D2h, D3h, D4h, D5h, D6h, Td, Oh, Ih,
    C_inf_v, D_inf_h
)

__version__ = "1.1.0"
__author__ = "J.R. Neilson"
__email__ = "james.neilson@colostate.edu"

__all__ = [
    "PointGroup",
    "PointGroupFactory", 
    "get_point_group",
    # Point groups
    "C1", "Ci", "Cs", "C2", "C3", "C2v", "C3v", "C4v", "C5v", "C6v",
    "D2", "D3", "D2h", "D3h", "D4h", "D5h", "D6h", "Td", "Oh", "Ih",
    "C_inf_v", "D_inf_h"
]