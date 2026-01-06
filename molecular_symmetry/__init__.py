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
    # Direct product functions
    calculate_direct_product,
    direct_product_label,
    # Point group classes
    C1, Ci, Cs, C2, C3, C4, C5, C6, C7, C8,
    C2v, C3v, C4v, C5v, C6v, C2h, C3h, C4h, C5h, C6h,
    D2, D3, D4, D5, D2h, D3h, D4h, D5h, D6h,
    D2d, D3d, D4d, D5d, D6d,
    S4, S6, S8,
    T, Td, Th, O, Oh, I, Ih,
    Cinfv, Dinfh,
    # Double groups for half-integer spin systems
    C2vStar, OhStar, TdStar, D4hStar, D2Star, D3Star, D6Star, OStar
)

__version__ = "1.1.0"
__author__ = "J.R. Neilson"
__email__ = "james.neilson@colostate.edu"

__all__ = [
    "PointGroup",
    "PointGroupFactory", 
    "get_point_group",
    # Direct product functions
    "calculate_direct_product",
    "direct_product_label",
    # Point groups
    "C1", "Ci", "Cs", "C2", "C3", "C4", "C5", "C6", "C7", "C8",
    "C2v", "C3v", "C4v", "C5v", "C6v", "C2h", "C3h", "C4h", "C5h", "C6h",
    "D2", "D3", "D4", "D5", "D2h", "D3h", "D4h", "D5h", "D6h",
    "D2d", "D3d", "D4d", "D5d", "D6d",
    "S4", "S6", "S8",
    "T", "Td", "Th", "O", "Oh", "I", "Ih",
    "Cinfv", "Dinfh",
    # Double groups
    "C2vStar", "OhStar", "TdStar", "D4hStar", "D2Star", "D3Star", "D6Star", "OStar"
]