"""
Molecular Symmetry Analysis Package
====================================
A comprehensive package for character tables and symmetry analysis
of molecules using Schoenflies notation.

Usage:
    from molecular_symmetry import PointGroupFactory
    
    pg = PointGroupFactory.create('Oh')
    irreps = pg.reduce_representation([6, 0, 0, 2, 2, 0, 0, 0, 4, 2])
"""

import numpy as np
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple


class PointGroup(ABC):
    """Abstract base class for point groups."""
    
    def __init__(self):
        self.name = ""
        self.classes = []
        self.class_sizes = []
        self.order = 0
        self.irreps = {}
    
    @abstractmethod
    def _initialize_character_table(self):
        """Initialize the character table for this point group."""
        pass
    
    def reduce_representation(self, reducible_rep: List[float]) -> Dict[str, int]:
        """
        Reduce a reducible representation to irreducible representations.
        
        Args:
            reducible_rep: List of characters for each symmetry class
            
        Returns:
            Dictionary mapping irrep names to their coefficients
        """
        if len(reducible_rep) != len(self.classes):
            raise ValueError(
                f"Reducible representation must have {len(self.classes)} characters "
                f"for point group {self.name}, got {len(reducible_rep)}"
            )
        
        coefficients = {}
        
        for irrep_name, irrep_chars in self.irreps.items():
            # Reduction formula: a_i = (1/h) * Σ(n_c * χ_reducible * χ_irrep)
            coeff = sum(
                class_size * char_red * char_irrep
                for class_size, char_red, char_irrep in 
                zip(self.class_sizes, reducible_rep, irrep_chars)
            ) / self.order
            
            coeff_int = int(round(coeff))
            if abs(coeff_int) > 0:
                coefficients[irrep_name] = coeff_int
        
        return coefficients
    
    def print_character_table(self):
        """Print the character table in a formatted way."""
        print(f"\n{'='*80}")
        print(f"Point Group: {self.name}")
        print(f"Order: {self.order}")
        print(f"{'='*80}")
        
        # Header - combine class size with class label for printing
        header = f"{'Irrep':<10}"
        for cls, size in zip(self.classes, self.class_sizes):
            if size > 1:
                header += f"{size}{cls:<8}"
            else:
                header += f"{cls:<10}"
        print(header)
        print("-" * 80)
        
        # Character table rows
        for irrep_name, chars in self.irreps.items():
            row = f"{irrep_name:<10}"
            for char in chars:
                if isinstance(char, complex):
                    row += f"{char:<10}"
                elif isinstance(char, float):
                    row += f"{char:<10.3f}"
                else:
                    row += f"{char:<10}"
            print(row)
    
    def get_symmetry_label(self, reducible_rep: List[float]) -> str:
        """Get a formatted symmetry label from a reducible representation."""
        irreps = self.reduce_representation(reducible_rep)
        if not irreps:
            return "0"
        
        terms = []
        for irrep, coeff in irreps.items():
            if coeff == 1:
                terms.append(irrep)
            else:
                terms.append(f"{coeff}{irrep}")
        
        return " ⊕ ".join(terms)
    
    def direct_product(self, irrep1: str, irrep2: str) -> Dict[str, int]:
        """
        Calculate the direct product of two irreducible representations.
        
        The direct product Γ1 ⊗ Γ2 gives the symmetry of the product of 
        functions belonging to irreps Γ1 and Γ2.
        
        Args:
            irrep1: Name of first irreducible representation
            irrep2: Name of second irreducible representation
            
        Returns:
            Dictionary mapping irrep names to their coefficients in the direct product
            
        Example:
            >>> oh = get_point_group('Oh')
            >>> oh.direct_product('T1u', 'T2g')
            {'T1g': 1, 'T2g': 1, 'Eg': 1, 'A2g': 1}
        """
        if irrep1 not in self.irreps:
            raise ValueError(f"Irrep '{irrep1}' not found in {self.name} character table")
        if irrep2 not in self.irreps:
            raise ValueError(f"Irrep '{irrep2}' not found in {self.name} character table")
        
        # Get characters for both irreps
        chars1 = self.irreps[irrep1]
        chars2 = self.irreps[irrep2]
        
        # Calculate direct product characters: χ(Γ1 ⊗ Γ2) = χ(Γ1) * χ(Γ2)
        product_chars = []
        for c1, c2 in zip(chars1, chars2):
            if isinstance(c1, complex) or isinstance(c2, complex):
                product_chars.append(complex(c1) * complex(c2))
            else:
                product_chars.append(c1 * c2)
        
        # Reduce the direct product to irreducible representations
        return self.reduce_representation(product_chars)
    
    def get_direct_product_table(self) -> Dict[Tuple[str, str], Dict[str, int]]:
        """
        Generate the complete direct product table for this point group.
        
        Returns:
            Dictionary mapping (irrep1, irrep2) tuples to their direct product decomposition
            
        Example:
            >>> oh = get_point_group('Oh')
            >>> table = oh.get_direct_product_table()
            >>> table[('T1u', 'T2g')]
            {'T1g': 1, 'T2g': 1, 'Eg': 1, 'A2g': 1}
        """
        table = {}
        irrep_names = list(self.irreps.keys())
        
        for i, irrep1 in enumerate(irrep_names):
            for j, irrep2 in enumerate(irrep_names):
                if j >= i:  # Only calculate upper triangle to avoid duplicates
                    table[(irrep1, irrep2)] = self.direct_product(irrep1, irrep2)
        
        return table
    
    def print_direct_product_table(self):
        """Print the direct product table in a formatted way."""
        print(f"\n{'='*80}")
        print(f"Direct Product Table for Point Group: {self.name}")
        print(f"{'='*80}")
        
        irrep_names = list(self.irreps.keys())
        
        # Print header
        header = f"{'×':<8}"
        for irrep in irrep_names:
            header += f"{irrep:<12}"
        print(header)
        print("-" * 80)
        
        # Print table rows
        for i, irrep1 in enumerate(irrep_names):
            row = f"{irrep1:<8}"
            for j, irrep2 in enumerate(irrep_names):
                if j >= i:
                    product = self.direct_product(irrep1, irrep2)
                    product_str = " ⊕ ".join([f"{coeff}{name}" if coeff > 1 else name 
                                             for name, coeff in product.items()])
                    row += f"{product_str:<12}"
                else:
                    row += f"{'·':<12}"  # Symmetric, so use dot for lower triangle
            print(row)


class C1(PointGroup):
    """C1 point group - no symmetry."""
    
    def __init__(self):
        super().__init__()
        self.name = "C1"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E']
        self.class_sizes = [1]
        self.order = 1
        self.irreps = {
            'A': [1]
        }


class Cs(PointGroup):
    """Cs point group - plane of symmetry."""
    
    def __init__(self):
        super().__init__()
        self.name = "Cs"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'σh']
        self.class_sizes = [1, 1]
        self.order = 2
        self.irreps = {
            "A'": [1, 1],
            'A"': [1, -1]
        }


class Ci(PointGroup):
    """Ci point group - inversion center."""
    
    def __init__(self):
        super().__init__()
        self.name = "Ci"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'i']
        self.class_sizes = [1, 1]
        self.order = 2
        self.irreps = {
            'Ag': [1, 1],
            'Au': [1, -1]
        }


class C2(PointGroup):
    """C2 point group - two-fold rotation axis."""
    
    def __init__(self):
        super().__init__()
        self.name = "C2"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C2']
        self.class_sizes = [1, 1]
        self.order = 2
        self.irreps = {
            'A': [1, 1],
            'B': [1, -1]
        }


class C3(PointGroup):
    """C3 point group - three-fold rotation axis."""
    
    def __init__(self):
        super().__init__()
        self.name = "C3"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C3', 'C3²']
        self.class_sizes = [1, 1, 1]
        self.order = 3
        
        ε = np.exp(2j * np.pi / 3)  # cube root of unity
        self.irreps = {
            'A': [1, 1, 1],
            'Ea': [1, ε, ε.conjugate()],
            'Eb': [1, ε.conjugate(), ε]
        }


class C2v(PointGroup):
    """C2v point group - two-fold axis with two vertical mirror planes."""
    
    def __init__(self):
        super().__init__()
        self.name = "C2v"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C2', 'σv(xz)', 'σv(yz)']
        self.class_sizes = [1, 1, 1, 1]
        self.order = 4
        self.irreps = {
            'A1': [1, 1, 1, 1],
            'A2': [1, 1, -1, -1],
            'B1': [1, -1, 1, -1],
            'B2': [1, -1, -1, 1]
        }


class C3v(PointGroup):
    """C3v point group - three-fold axis with three vertical mirror planes."""
    
    def __init__(self):
        super().__init__()
        self.name = "C3v"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C3', 'σv']
        self.class_sizes = [1, 2, 3]
        self.order = 6
        self.irreps = {
            'A1': [1, 1, 1],
            'A2': [1, 1, -1],
            'E': [2, -1, 0]
        }


class C4v(PointGroup):
    """C4v point group - four-fold axis with vertical mirror planes."""
    
    def __init__(self):
        super().__init__()
        self.name = "C4v"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C4', 'C2', 'σv', 'σd']
        self.class_sizes = [1, 2, 1, 2, 2]
        self.order = 8
        self.irreps = {
            'A1': [1, 1, 1, 1, 1],
            'A2': [1, 1, 1, -1, -1],
            'B1': [1, -1, 1, 1, -1],
            'B2': [1, -1, 1, -1, 1],
            'E': [2, 0, -2, 0, 0]
        }


class C5v(PointGroup):
    """C5v point group - five-fold axis with vertical mirror planes."""
    
    def __init__(self):
        super().__init__()
        self.name = "C5v"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C5', 'C5²', 'σv']
        self.class_sizes = [1, 2, 2, 5]
        self.order = 10
        
        φ = (1 + np.sqrt(5)) / 2  # golden ratio
        self.irreps = {
            'A1': [1, 1, 1, 1],
            'A2': [1, 1, 1, -1],
            'E1': [2, 2*np.cos(2*np.pi/5), 2*np.cos(4*np.pi/5), 0],
            'E2': [2, 2*np.cos(4*np.pi/5), 2*np.cos(2*np.pi/5), 0]
        }


class C6v(PointGroup):
    """C6v point group - six-fold axis with vertical mirror planes."""
    
    def __init__(self):
        super().__init__()
        self.name = "C6v"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C6', 'C3', 'C2', 'σv', 'σd']
        self.class_sizes = [1, 2, 2, 1, 3, 3]
        self.order = 12
        self.irreps = {
            'A1': [1, 1, 1, 1, 1, 1],
            'A2': [1, 1, 1, 1, -1, -1],
            'B1': [1, -1, 1, -1, 1, -1],
            'B2': [1, -1, 1, -1, -1, 1],
            'E1': [2, 1, -1, -2, 0, 0],
            'E2': [2, -1, -1, 2, 0, 0]
        }


class D2(PointGroup):
    """D2 point group - three perpendicular two-fold axes."""
    
    def __init__(self):
        super().__init__()
        self.name = "D2"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C2(z)', 'C2(y)', 'C2(x)']
        self.class_sizes = [1, 1, 1, 1]
        self.order = 4
        self.irreps = {
            'A': [1, 1, 1, 1],
            'B1': [1, 1, -1, -1],
            'B2': [1, -1, 1, -1],
            'B3': [1, -1, -1, 1]
        }


class D2h(PointGroup):
    """D2h point group - D2 with inversion center."""
    
    def __init__(self):
        super().__init__()
        self.name = "D2h"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C2(z)', 'C2(y)', 'C2(x)', 'i', 'σ(xy)', 'σ(xz)', 'σ(yz)']
        self.class_sizes = [1, 1, 1, 1, 1, 1, 1, 1]
        self.order = 8
        self.irreps = {
            'Ag': [1, 1, 1, 1, 1, 1, 1, 1],
            'B1g': [1, 1, -1, -1, 1, 1, -1, -1],
            'B2g': [1, -1, 1, -1, 1, -1, 1, -1],
            'B3g': [1, -1, -1, 1, 1, -1, -1, 1],
            'Au': [1, 1, 1, 1, -1, -1, -1, -1],
            'B1u': [1, 1, -1, -1, -1, -1, 1, 1],
            'B2u': [1, -1, 1, -1, -1, 1, -1, 1],
            'B3u': [1, -1, -1, 1, -1, 1, 1, -1]
        }


class D3(PointGroup):
    """D3 point group - three-fold axis with perpendicular two-fold axes."""
    
    def __init__(self):
        super().__init__()
        self.name = "D3"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C3', 'C2']
        self.class_sizes = [1, 2, 3]
        self.order = 6
        self.irreps = {
            'A1': [1, 1, 1],
            'A2': [1, 1, -1],
            'E': [2, -1, 0]
        }


class D3h(PointGroup):
    """D3h point group - D3 with horizontal mirror plane."""
    
    def __init__(self):
        super().__init__()
        self.name = "D3h"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C3', 'C2', 'σh', 'S3', 'σv']
        self.class_sizes = [1, 2, 3, 1, 2, 3]
        self.order = 12
        self.irreps = {
            "A1'": [1, 1, 1, 1, 1, 1],
            "A2'": [1, 1, -1, 1, 1, -1],
            "E'": [2, -1, 0, 2, -1, 0],
            'A1"': [1, 1, 1, -1, -1, -1],
            'A2"': [1, 1, -1, -1, -1, 1],
            'E"': [2, -1, 0, -2, 1, 0]
        }


class D4h(PointGroup):
    """D4h point group - four-fold axis with horizontal mirror plane."""
    
    def __init__(self):
        super().__init__()
        self.name = "D4h"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C4', 'C2', 'C2\'', 'C2"', 'i', 'S4', 'σh', 'σv', 'σd']
        self.class_sizes = [1, 2, 1, 2, 2, 1, 2, 1, 2, 2]
        self.order = 16
        self.irreps = {
            'A1g': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            'A2g': [1, 1, 1, -1, -1, 1, 1, 1, -1, -1],
            'B1g': [1, -1, 1, 1, -1, 1, -1, 1, 1, -1],
            'B2g': [1, -1, 1, -1, 1, 1, -1, 1, -1, 1],
            'Eg': [2, 0, -2, 0, 0, 2, 0, -2, 0, 0],
            'A1u': [1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
            'A2u': [1, 1, 1, -1, -1, -1, -1, -1, 1, 1],
            'B1u': [1, -1, 1, 1, -1, -1, 1, -1, -1, 1],
            'B2u': [1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
            'Eu': [2, 0, -2, 0, 0, -2, 0, 2, 0, 0]
        }


class D5h(PointGroup):
    """D5h point group - five-fold axis with horizontal mirror plane."""
    
    def __init__(self):
        super().__init__()
        self.name = "D5h"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C5', 'C5²', 'C2', 'σh', 'S5', 'S5³', 'σv']
        self.class_sizes = [1, 2, 2, 5, 1, 2, 2, 5]
        self.order = 20
        
        c1 = 2*np.cos(2*np.pi/5)
        c2 = 2*np.cos(4*np.pi/5)
        
        self.irreps = {
            "A1'": [1, 1, 1, 1, 1, 1, 1, 1],
            "A2'": [1, 1, 1, -1, 1, 1, 1, -1],
            "E1'": [2, c1, c2, 0, 2, c1, c2, 0],
            "E2'": [2, c2, c1, 0, 2, c2, c1, 0],
            'A1"': [1, 1, 1, 1, -1, -1, -1, -1],
            'A2"': [1, 1, 1, -1, -1, -1, -1, 1],
            'E1"': [2, c1, c2, 0, -2, -c1, -c2, 0],
            'E2"': [2, c2, c1, 0, -2, -c2, -c1, 0]
        }


class D6h(PointGroup):
    """D6h point group - six-fold axis with horizontal mirror plane."""
    
    def __init__(self):
        super().__init__()
        self.name = "D6h"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C6', 'C3', 'C2', 'C2\'', 'C2"', 
                       'i', 'S3', 'S6', 'σh', 'σd', 'σv']
        self.class_sizes = [1, 2, 2, 1, 3, 3, 1, 2, 2, 1, 3, 3]
        self.order = 24
        self.irreps = {
            'A1g': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            'A2g': [1, 1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1],
            'B1g': [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
            'B2g': [1, -1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1],
            'E1g': [2, 1, -1, -2, 0, 0, 2, 1, -1, -2, 0, 0],
            'E2g': [2, -1, -1, 2, 0, 0, 2, -1, -1, 2, 0, 0],
            'A1u': [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1],
            'A2u': [1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1],
            'B1u': [1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1],
            'B2u': [1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1],
            'E1u': [2, 1, -1, -2, 0, 0, -2, -1, 1, 2, 0, 0],
            'E2u': [2, -1, -1, 2, 0, 0, -2, 1, 1, -2, 0, 0]
        }


class Td(PointGroup):
    """Td point group - tetrahedral symmetry."""
    
    def __init__(self):
        super().__init__()
        self.name = "Td"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C3', 'C2', 'S4', 'σd']
        self.class_sizes = [1, 8, 3, 6, 6]
        self.order = 24
        self.irreps = {
            'A1': [1, 1, 1, 1, 1],
            'A2': [1, 1, 1, -1, -1],
            'E': [2, -1, 2, 0, 0],
            'T1': [3, 0, -1, 1, -1],
            'T2': [3, 0, -1, -1, 1]
        }


class Oh(PointGroup):
    """Oh point group - octahedral symmetry."""
    
    def __init__(self):
        super().__init__()
        self.name = "Oh"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C3', 'C2', 'C4', 'C2\'', 'i', 'S4', 'S6', 'σh', 'σd']
        self.class_sizes = [1, 8, 6, 6, 3, 1, 6, 8, 3, 6]
        self.order = 48
        self.irreps = {
            'A1g': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            'A2g': [1, 1, -1, -1, 1, 1, -1, 1, 1, -1],
            'Eg': [2, -1, 0, 0, 2, 2, 0, -1, 2, 0],
            'T1g': [3, 0, -1, 1, -1, 3, 1, 0, -1, -1],
            'T2g': [3, 0, 1, -1, -1, 3, -1, 0, -1, 1],
            'A1u': [1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
            'A2u': [1, 1, -1, -1, 1, -1, 1, -1, -1, 1],
            'Eu': [2, -1, 0, 0, 2, -2, 0, 1, -2, 0],
            'T1u': [3, 0, -1, 1, -1, -3, -1, 0, 1, 1],
            'T2u': [3, 0, 1, -1, -1, -3, 1, 0, 1, -1]
        }


class Ih(PointGroup):
    """Ih point group - icosahedral symmetry."""
    
    def __init__(self):
        super().__init__()
        self.name = "Ih"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C5', 'C5²', 'C3', 'C2', 
                       'i', 'S10', 'S10³', 'S6', 'σ']
        self.class_sizes = [1, 12, 12, 20, 15, 1, 12, 12, 20, 15]
        self.order = 120
        
        τ = (1 + np.sqrt(5)) / 2  # golden ratio
        
        self.irreps = {
            'Ag': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            'T1g': [3, τ, -1/τ, 0, -1, 3, -1/τ, τ, 0, -1],
            'T2g': [3, -1/τ, τ, 0, -1, 3, τ, -1/τ, 0, -1],
            'Gg': [4, -1, -1, 1, 0, 4, -1, -1, 1, 0],
            'Hg': [5, 0, 0, -1, 1, 5, 0, 0, -1, 1],
            'Au': [1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
            'T1u': [3, τ, -1/τ, 0, -1, -3, 1/τ, -τ, 0, 1],
            'T2u': [3, -1/τ, τ, 0, -1, -3, -τ, 1/τ, 0, 1],
            'Gu': [4, -1, -1, 1, 0, -4, 1, 1, -1, 0],
            'Hu': [5, 0, 0, -1, 1, -5, 0, 0, 1, -1]
        }


class Cinfv(PointGroup):
    """C∞v point group - linear molecules (heteronuclear)."""
    
    def __init__(self):
        super().__init__()
        self.name = "C∞v"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        # Simplified representation for linear molecules
        self.classes = ['E', 'C∞φ', 'σv']
        self.class_sizes = [1, 1, 1]  # Conceptual
        self.order = float('inf')
        self.irreps = {
            'Σ+': [1, 1, 1],
            'Σ-': [1, 1, -1],
            'Π': [2, 2, 0],
            'Δ': [2, 2, 0],
            'Φ': [2, 2, 0]
        }
    
    def reduce_representation(self, reducible_rep: List[float]) -> Dict[str, int]:
        """Linear molecules require special handling."""
        print("Note: C∞v is an infinite group. Reduction is approximate.")
        return super().reduce_representation(reducible_rep)


class Dinfh(PointGroup):
    """D∞h point group - linear molecules (homonuclear)."""
    
    def __init__(self):
        super().__init__()
        self.name = "D∞h"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        # Simplified representation for linear molecules
        self.classes = ['E', 'C∞φ', 'σv', 'i', 'S∞φ', 'C2']
        self.class_sizes = [1, 1, 1, 1, 1, 1]  # Conceptual
        self.order = float('inf')
        self.irreps = {
            'Σg+': [1, 1, 1, 1, 1, 1],
            'Σg-': [1, 1, -1, 1, 1, -1],
            'Πg': [2, 2, 0, 2, 2, 0],
            'Δg': [2, 2, 0, 2, 2, 0],
            'Σu+': [1, 1, 1, -1, -1, -1],
            'Σu-': [1, 1, -1, -1, -1, 1],
            'Πu': [2, 2, 0, -2, -2, 0],
            'Δu': [2, 2, 0, -2, -2, 0]
        }
    
    def reduce_representation(self, reducible_rep: List[float]) -> Dict[str, int]:
        """Linear molecules require special handling."""
        print("Note: D∞h is an infinite group. Reduction is approximate.")
        return super().reduce_representation(reducible_rep)


class PointGroupFactory:
    """Factory class to create point group objects."""
    
    _point_groups = {
        'C1': C1,
        'Ci': Ci,
        'Cs': Cs,
        'C2': C2,
        'C3': C3,
        'C2v': C2v,
        'C3v': C3v,
        'C4v': C4v,
        'C5v': C5v,
        'C6v': C6v,
        'D2': D2,
        'D2h': D2h,
        'D3': D3,
        'D3h': D3h,
        'D4h': D4h,
        'D5h': D5h,
        'D6h': D6h,
        'Td': Td,
        'Oh': Oh,
        'Ih': Ih,
        'Cinfv': Cinfv,
        'C∞v': Cinfv,
        'Dinfh': Dinfh,
        'D∞h': Dinfh
    }
    
    @classmethod
    def create(cls, point_group_name: str) -> PointGroup:
        """
        Create a point group object by name.
        
        Args:
            point_group_name: Name of the point group (e.g., 'Oh', 'Td', 'D4h')
            
        Returns:
            PointGroup object
            
        Raises:
            ValueError: If point group name is not recognized
        """
        if point_group_name not in cls._point_groups:
            available = ', '.join(sorted(cls._point_groups.keys()))
            raise ValueError(
                f"Unknown point group '{point_group_name}'. "
                f"Available point groups: {available}"
            )
        
        return cls._point_groups[point_group_name]()
    
    @classmethod
    def list_available(cls) -> List[str]:
        """Return a list of all available point groups."""
        return sorted(cls._point_groups.keys())


# Convenience function
def get_point_group(name: str) -> PointGroup:
    """
    Convenience function to get a point group.
    
    Args:
        name: Point group name (e.g., 'Oh', 'Td', 'C2v')
        
    Returns:
        PointGroup object
    """
    return PointGroupFactory.create(name)


def calculate_direct_product(point_group_name: str, irrep1: str, irrep2: str) -> Dict[str, int]:
    """
    Convenience function to calculate direct product of two irreps.
    
    Args:
        point_group_name: Name of the point group (e.g., 'Oh', 'Td')
        irrep1: First irreducible representation
        irrep2: Second irreducible representation
        
    Returns:
        Dictionary mapping irrep names to coefficients in the direct product
        
    Example:
        >>> calculate_direct_product('Oh', 'T1u', 'T2g')
        {'T1g': 1, 'T2g': 1, 'Eg': 1, 'A2g': 1}
    """
    pg = get_point_group(point_group_name)
    return pg.direct_product(irrep1, irrep2)


def direct_product_label(point_group_name: str, irrep1: str, irrep2: str) -> str:
    """
    Get formatted label for direct product of two irreps.
    
    Args:
        point_group_name: Name of the point group
        irrep1: First irreducible representation
        irrep2: Second irreducible representation
        
    Returns:
        Formatted string showing the direct product decomposition
        
    Example:
        >>> direct_product_label('Oh', 'T1u', 'T2g')
        'T1g ⊕ T2g ⊕ Eg ⊕ A2g'
    """
    product = calculate_direct_product(point_group_name, irrep1, irrep2)
    terms = []
    for irrep, coeff in product.items():
        if coeff == 1:
            terms.append(irrep)
        else:
            terms.append(f"{coeff}{irrep}")
    return " ⊕ ".join(terms)


if __name__ == "__main__":
    # Demo of the package
    print("Molecular Symmetry Analysis Package")
    print("=" * 80)
    print("\nAvailable point groups:")
    print(", ".join(PointGroupFactory.list_available()))
    
    # Example with Oh
    print("\n\nExample: Octahedral ML6 Complex")
    oh = get_point_group('Oh')
    oh.print_character_table()
    
    print("\n\nSigma bonding orbitals (6 ligands):")
    sigma_rep = [6, 0, 0, 2, 2, 0, 0, 0, 4, 2]
    print(f"Γ_σ = {sigma_rep}")
    print(f"Reduction: {oh.get_symmetry_label(sigma_rep)}")
    
    print("\nMetal d orbitals:")
    d_rep = [5, 2, 1, 1, 1, 5, 1, 2, 1, 1]
    print(f"Γ_d = {d_rep}")
    print(f"Reduction: {oh.get_symmetry_label(d_rep)}")
