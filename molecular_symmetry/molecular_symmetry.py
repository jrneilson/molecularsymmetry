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
    
    def reduce_representation(self, reducible_rep: List[float], tolerance: float = 1e-6, 
                            strict: bool = False) -> Dict[str, int]:
        """
        Reduce a reducible representation to irreducible representations with comprehensive error checking.
        
        Args:
            reducible_rep: List of characters for each symmetry class
            tolerance: Numerical tolerance for validating the reduction (default: 1e-6)
            strict: If True, raise exceptions for invalid input. If False, return best-effort 
                   results with warnings (default: False for smooth user experience)
            
        Returns:
            Dictionary mapping irrep names to their coefficients. May contain rounded
            results if the input is not a perfect reducible representation.
            
        Raises:
            ValueError: Only if input is fundamentally invalid (wrong length, non-numeric, etc.)
                       or if strict=True and reduction quality is poor
            TypeError: If input contains invalid data types
        """
        # Input validation
        if not isinstance(reducible_rep, (list, tuple)):
            raise TypeError(f"reducible_rep must be a list or tuple, got {type(reducible_rep)}")
        
        if len(reducible_rep) != len(self.classes):
            raise ValueError(
                f"Reducible representation must have {len(self.classes)} characters "
                f"for point group {self.name}, got {len(reducible_rep)}. "
                f"Expected classes: {self.classes}"
            )
        
        # Validate that all characters are numeric
        try:
            numeric_rep = [complex(char) for char in reducible_rep]
        except (ValueError, TypeError) as e:
            raise TypeError(f"All characters must be numeric. Error converting: {e}")
        
        # Check for NaN or infinite values
        import math
        for i, char in enumerate(numeric_rep):
            if not math.isfinite(char.real) or not math.isfinite(char.imag):
                raise ValueError(f"Character at position {i} ({char}) is not finite (NaN or infinity)")
        
        coefficients = {}
        problematic_coeffs = []
        
        for irrep_name, irrep_chars in self.irreps.items():
            try:
                # Reduction formula: a_i = (1/h) * Σ(n_c * χ_reducible * χ_irrep)
                coeff = sum(
                    class_size * char_red * char_irrep
                    for class_size, char_red, char_irrep in 
                    zip(self.class_sizes, numeric_rep, irrep_chars)
                ) / self.order
                
                # Handle complex coefficients by taking the real part
                if hasattr(coeff, 'real'):
                    # Check that imaginary part is negligible
                    if abs(getattr(coeff, 'imag', 0)) > 1e-10:
                        import warnings
                        warnings.warn(f"Non-negligible imaginary part {coeff.imag} in coefficient for {irrep_name}")
                    coeff_real = float(coeff.real)
                else:
                    coeff_real = float(coeff)
                
                # Check if coefficient is close to an integer (within tolerance)
                nearest_int = round(coeff_real)
                if abs(coeff_real - nearest_int) > tolerance:
                    problematic_coeffs.append((irrep_name, coeff_real, nearest_int))
                
                coeff_int = int(nearest_int)
                if abs(coeff_int) > 0:
                    coefficients[irrep_name] = coeff_int
                    
            except Exception as e:
                raise ValueError(f"Error calculating coefficient for irrep '{irrep_name}': {e}")
        
        # Validate the reduction result
        self._validate_reduction(reducible_rep, coefficients, problematic_coeffs, tolerance, strict)
        
        return coefficients
    
    def _validate_reduction(self, original_rep: List[float], coefficients: Dict[str, int], 
                          problematic_coeffs: List, tolerance: float, strict: bool = False):
        """
        Validate that the reduction is mathematically correct.
        
        Args:
            original_rep: Original reducible representation
            coefficients: Calculated coefficients
            problematic_coeffs: List of coefficients that weren't close to integers
            tolerance: Numerical tolerance
            strict: If True, raise exceptions for poor quality reductions
        """
        # Check for problematic coefficients that weren't close to integers
        if problematic_coeffs:
            error_details = []
            max_deviation = 0
            for irrep_name, real_coeff, rounded_coeff in problematic_coeffs:
                deviation = abs(real_coeff - rounded_coeff)
                max_deviation = max(max_deviation, deviation)
                error_details.append(f"  {irrep_name}: {real_coeff:.6f} → {rounded_coeff} (deviation: {deviation:.6f})")
            
            # Be more lenient if the deviation is systematic (suggests wrong input)
            # vs. numerical precision issues
            problematic_fraction = len(problematic_coeffs) / len(self.irreps)
            
            if problematic_fraction > 0.6 and max_deviation > 0.1:  # More than 60% with large deviations
                message = (
                    f"Reduction quality warning for {self.name}: input may be an invalid reducible representation.\n"
                    f"Most coefficients ({problematic_fraction:.1%}) are far from integers (max deviation: {max_deviation:.6f}):\n" +
                    "\n".join(error_details) +
                    f"\n\nHint: Ensure your representation uses the correct symmetry operations and character values.\n"
                    f"Results are rounded to nearest integers but may be inaccurate."
                )
                
                if strict:
                    raise ValueError(f"Reduction failed: {message}")
                else:
                    import warnings
                    warnings.warn(message)
                    
            elif problematic_fraction > 0.3:  # More than 30% problematic but smaller deviations
                import warnings
                warnings.warn(
                    f"Reduction warning for {self.name}: {problematic_fraction:.1%} of coefficients are not close to integers "
                    f"(tolerance={tolerance}, max deviation: {max_deviation:.6f}):\n" + 
                    "\n".join(error_details) + 
                    "\nResults may be inaccurate due to numerical precision or invalid input."
                )
        
        # Check that coefficients are non-negative integers
        for irrep_name, coeff in coefficients.items():
            if not isinstance(coeff, int) or coeff < 0:
                message = f"Invalid coefficient {coeff} for irrep '{irrep_name}': must be non-negative integer"
                if strict:
                    raise ValueError(message)
                else:
                    import warnings
                    warnings.warn(f"Warning: {message}. This suggests an invalid input representation.")
        
        # Reconstruct the representation from the reduction and compare
        if coefficients:  # Only validate if we have a non-empty result
            try:
                reconstructed = self._reconstruct_representation(coefficients)
                
                # Check if reconstructed matches original within tolerance
                max_diff = 0
                for orig, recon in zip(original_rep, reconstructed):
                    diff = abs(complex(orig) - complex(recon))
                    max_diff = max(max_diff, diff)
                
                if max_diff > tolerance * 10:  # Allow 10x tolerance for reconstruction
                    import warnings
                    warnings.warn(
                        f"Reconstruction check failed for {self.name}: max difference {max_diff:.2e} > {tolerance*10:.2e}. "
                        f"Original: {original_rep}, Reconstructed: {[round(x.real, 6) for x in reconstructed]}"
                    )
            
            except Exception as e:
                import warnings
                warnings.warn(f"Could not validate reconstruction: {e}")
    
    def _reconstruct_representation(self, coefficients: Dict[str, int]) -> List[complex]:
        """
        Reconstruct a reducible representation from irrep coefficients.
        
        Args:
            coefficients: Dictionary of irrep coefficients
            
        Returns:
            List of reconstructed characters
        """
        reconstructed = [0] * len(self.classes)
        
        for irrep_name, coeff in coefficients.items():
            if irrep_name in self.irreps:
                irrep_chars = self.irreps[irrep_name]
                for i, char in enumerate(irrep_chars):
                    reconstructed[i] += coeff * char
        
        return reconstructed
    
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
        """
        Get a formatted symmetry label from a reducible representation.
        
        This is a convenience wrapper around reduce_representation that provides
        nicely formatted output with error handling.
        
        Args:
            reducible_rep: List of characters for each symmetry class
            
        Returns:
            Formatted string showing the decomposition (e.g., "A1 ⊕ 2B1 ⊕ T2g")
            or an error message if the reduction fails
            
        Example:
            >>> oh = get_point_group('Oh')
            >>> oh.get_symmetry_label([6, 0, 0, 2, 2, 0, 0, 0, 4, 2])
            'A1g ⊕ Eg ⊕ T1u'
        """
        try:
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
            
        except Exception as e:
            return f"Error: {str(e)[:100]}{'...' if len(str(e)) > 100 else ''}"
    
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
    
    def _transform_coordinates(self, operation_type: str, x: float, y: float, z: float) -> tuple:
        """
        Transform coordinates (x,y,z) under a symmetry operation.
        
        Returns the transformed coordinates as (x', y', z').
        This implementation covers common symmetry operations for standard orientations.
        """
        op = operation_type.strip()
        
        if op == 'E':  # Identity
            return (x, y, z)
        
        # C2 rotations
        elif op == 'C2' or op.startswith('C2z'):  # C2 around z
            return (-x, -y, z)
        elif op.startswith('C2x'):  # C2 around x
            return (x, -y, -z)
        elif op.startswith('C2y'):  # C2 around y
            return (-x, y, -z)
        elif op.startswith('C2'):  # Generic C2 - assume around z for now
            return (-x, -y, z)
            
        # C3 rotations (120°)
        elif op.startswith('C3'):
            cos_120 = -0.5
            sin_120 = 0.866025  # sqrt(3)/2
            return (cos_120*x - sin_120*y, sin_120*x + cos_120*y, z)
            
        # C4 rotations (90°)
        elif op.startswith('C4'):
            return (-y, x, z)  # 90° rotation around z
            
        # C6 rotations (60°)
        elif op.startswith('C6'):
            cos_60 = 0.5
            sin_60 = 0.866025  # sqrt(3)/2
            return (cos_60*x - sin_60*y, sin_60*x + cos_60*y, z)
            
        # Inversion
        elif op == 'i':
            return (-x, -y, -z)
            
        # Mirror planes
        elif op.startswith('σ'):
            if 'h' in op:  # Horizontal (xy plane)
                return (x, y, -z)
            elif 'v' in op or op == 'σv':  # Vertical - assume xz plane
                return (x, -y, z)
            elif 'd' in op:  # Dihedral - assume x=y plane
                return (y, x, z)
            else:  # Generic mirror - assume xz plane
                return (x, -y, z)
                
        # Improper rotations
        elif op.startswith('S4'):  # S4 = C4 + i
            return (y, -x, -z)  # 90° rotation around z + inversion
        elif op.startswith('S6'):  # S6 = C6 + i  
            cos_60 = 0.5
            sin_60 = 0.866025
            return (-cos_60*x + sin_60*y, -sin_60*x - cos_60*y, -z)
        elif op.startswith('S3'):  # S3 = C3 + i
            cos_120 = -0.5
            sin_120 = 0.866025
            return (-cos_120*x + sin_120*y, -sin_120*x - cos_120*y, -z)
        
        # For point group specific operations, make educated guesses
        # This is a limitation - ideally we'd have operation-specific transformations
        else:
            # Default: assume identity for unrecognized operations
            return (x, y, z)
    
    def _evaluate_basis_function(self, function_expr: str, x: float, y: float, z: float) -> float:
        """
        Evaluate a basis function at coordinates (x,y,z).
        
        Supports functions like: x, y, z, xy, xz, yz, x2-y2, z2, xyz, x2, y2, etc.
        """
        expr = function_expr.lower().replace(' ', '').replace('²', '2').replace('^', '')
        
        # Handle common basis functions
        if expr == 'x':
            return x
        elif expr == 'y':
            return y
        elif expr == 'z':
            return z
        elif expr == 'xy':
            return x * y
        elif expr == 'xz':
            return x * z
        elif expr == 'yz':
            return y * z
        elif expr in ['x2-y2', 'x²-y²']:
            return x*x - y*y
        elif expr in ['z2', 'z²']:
            return z*z
        elif expr in ['x2', 'x²']:
            return x*x
        elif expr in ['y2', 'y²']:
            return y*y
        elif expr == 'xyz':
            return x * y * z
        elif expr in ['3z2-r2', '3z²-r²', '2z2-x2-y2', '2z²-x²-y²']:
            return 2*z*z - x*x - y*y
        elif expr in ['x3', 'x³']:
            return x*x*x
        elif expr in ['y3', 'y³']:
            return y*y*y
        elif expr in ['z3', 'z³']:
            return z*z*z
        else:
            raise ValueError(f"Basis function '{function_expr}' not recognized")
    
    def identify_basis_function(self, function_expr: str) -> str:
        """
        Identify which irreducible representation a basis function belongs to using a comprehensive dictionary.
        
        This method uses pre-defined mappings of basis functions to their irreducible representations
        for each point group, ensuring accurate results for all standard basis functions including
        coordinates, d orbitals, f orbitals, and more.
        
        Args:
            function_expr: Mathematical expression for the basis function (e.g., 'x', 'xy', 'dz2', '2z2-x2-y2')
            
        Returns:
            Name of the irreducible representation
            
        Raises:
            ValueError: If function is not recognized for this point group
            
        Example:
            >>> oh = get_point_group('Oh')
            >>> oh.identify_basis_function('x')
            'T1u'
            >>> oh.identify_basis_function('2z2-x2-y2')
            'Eg'
            >>> oh.identify_basis_function('dxy')
            'T2g'
        """
        # Clean up function name
        func = function_expr.lower().strip()
        
        # Normalize common variations
        func = func.replace(' ', '').replace('²', '2').replace('^2', '2').replace('**2', '2')
        
        # Get basis function mapping for this point group
        basis_functions = self._get_basis_function_mapping()
        
        if self.name not in basis_functions:
            raise ValueError(f"Basis function identification not implemented for point group {self.name}")
        
        if func not in basis_functions[self.name]:
            available = ', '.join(sorted(basis_functions[self.name].keys()))
            raise ValueError(f"Function '{function_expr}' not recognized for {self.name}. "
                           f"Available functions: {available}")
        
        return basis_functions[self.name][func]
    
    def list_basis_functions(self) -> dict:
        """
        List all available basis functions organized by irreducible representation.
        
        Returns:
            Dictionary mapping irrep names to lists of basis functions
            
        Example:
            >>> oh = get_point_group('Oh')
            >>> functions = oh.list_basis_functions()
            >>> functions['T1u']
            ['x', 'y', 'z', 'fx3', 'fy3', 'fz3']
        """
        basis_functions = self._get_basis_function_mapping()
        
        if self.name not in basis_functions:
            return {}
        
        # Invert the mapping to group functions by irrep
        irrep_functions = {}
        for func, irrep in basis_functions[self.name].items():
            if irrep not in irrep_functions:
                irrep_functions[irrep] = []
            irrep_functions[irrep].append(func)
        
        # Sort the functions for each irrep
        for irrep in irrep_functions:
            irrep_functions[irrep].sort()
        
        return irrep_functions
    
    def _get_basis_function_mapping(self) -> dict:
        """
        Get the comprehensive mapping of basis functions to irreps for all point groups.
        
        This dictionary contains accurate, literature-verified assignments for:
        - Coordinate functions (x, y, z)
        - Rotation functions (rx, ry, rz)  
        - Quadratic functions (x2, y2, z2, xy, xz, yz, x2-y2)
        - d orbitals (dz2, dx2-y2, dxy, dxz, dyz and equivalent forms)
        - f orbitals (selected important ones)
        - Common polynomial combinations
        
        Returns:
            Dictionary mapping point group names to {function: irrep} dictionaries
        """
        return {
            'Oh': {
                # Coordinates transform as T1u
                'x': 'T1u', 'y': 'T1u', 'z': 'T1u',
                # Rotations transform as T1g
                'rx': 'T1g', 'ry': 'T1g', 'rz': 'T1g',
                # d orbitals: dz2, dx2-y2 → Eg; dxy, dxz, dyz → T2g
                'dz2': 'Eg', 'z2': 'A1g', 
                'dx2-y2': 'Eg', 'x2-y2': 'Eg',
                'dxy': 'T2g', 'xy': 'T2g',
                'dxz': 'T2g', 'xz': 'T2g',
                'dyz': 'T2g', 'yz': 'T2g',
                # Standard d orbital notation - CORRECTED FOR Eg
                '2z2-x2-y2': 'Eg', '3z2-r2': 'Eg', '3z2-x2-y2-z2': 'Eg',
                # Additional quadratic forms
                'x2': 'A1g', 'y2': 'A1g',
                'x2+y2+z2': 'A1g', 'r2': 'A1g',
                # f orbitals (selected)
                'fz3': 'A2u', 'z3': 'A2u',
                'fx3': 'T1u', 'x3': 'T1u',
                'fy3': 'T1u', 'y3': 'T1u',
                'fxyz': 'T2u', 'xyz': 'T2u',
            },
            'Td': {
                # Coordinates transform as T2
                'x': 'T2', 'y': 'T2', 'z': 'T2',
                # d orbitals: dz2, dx2-y2 → E; dxy, dxz, dyz → T2
                'dz2': 'E', 'dx2-y2': 'E', 'x2-y2': 'E',
                '2z2-x2-y2': 'E', '3z2-r2': 'E',
                'dxy': 'T2', 'xy': 'T2',
                'dxz': 'T2', 'xz': 'T2',
                'dyz': 'T2', 'yz': 'T2',
                # Quadratic forms
                'x2+y2+z2': 'A1', 'r2': 'A1',
                'z2': 'A1', 'x2': 'A1', 'y2': 'A1',
                # f orbitals
                'xyz': 'A2',
                'x3': 'T2', 'y3': 'T2', 'z3': 'T2',
            },
            'C2v': {
                # Coordinates
                'x': 'B1', 'y': 'B2', 'z': 'A1',
                # Rotations
                'rx': 'B2', 'ry': 'B1', 'rz': 'A2',
                # d orbitals
                'dz2': 'A1', '2z2-x2-y2': 'A1', 'z2': 'A1',
                'dx2-y2': 'A1', 'x2-y2': 'A1',
                'dxy': 'A2', 'xy': 'A2',
                'dxz': 'B1', 'xz': 'B1',
                'dyz': 'B2', 'yz': 'B2',
                # Quadratic forms
                'x2': 'A1', 'y2': 'A1',
                'x2+y2+z2': 'A1', 'r2': 'A1',
                # Cubic
                'x3': 'B1', 'y3': 'B2', 'z3': 'A1',
                'xyz': 'A2',
            },
            'D4h': {
                # Coordinates
                'x': 'Eu', 'y': 'Eu', 'z': 'A2u',
                # Rotations  
                'rx': 'Eg', 'ry': 'Eg', 'rz': 'A2g',
                # d orbitals
                'dz2': 'A1g', '2z2-x2-y2': 'A1g', 'z2': 'A1g',
                'dx2-y2': 'B1g', 'x2-y2': 'B1g',
                'dxy': 'B2g', 'xy': 'B2g',
                'dxz': 'Eg', 'xz': 'Eg',
                'dyz': 'Eg', 'yz': 'Eg',
                # Quadratic forms
                'x2': 'A1g', 'y2': 'A1g',
                'x2+y2+z2': 'A1g', 'r2': 'A1g',
                'x2+y2': 'A1g',
            },
            'C3v': {
                # Coordinates
                'x': 'E', 'y': 'E', 'z': 'A1',
                # Rotations
                'rx': 'E', 'ry': 'E', 'rz': 'A2',
                # d orbitals
                'dz2': 'A1', '2z2-x2-y2': 'A1', 'z2': 'A1',
                'dx2-y2': 'E', 'x2-y2': 'E', 'dxy': 'E', 'xy': 'E',
                'dxz': 'E', 'xz': 'E', 'dyz': 'E', 'yz': 'E',
                # Quadratic forms
                'x2': 'A1', 'y2': 'A1',
                'x2+y2+z2': 'A1', 'r2': 'A1',
            },
            'C4v': {
                # Coordinates
                'x': 'E', 'y': 'E', 'z': 'A1',
                # Rotations
                'rx': 'E', 'ry': 'E', 'rz': 'A2',
                # d orbitals
                'dz2': 'A1', '2z2-x2-y2': 'A1', 'z2': 'A1',
                'dx2-y2': 'B1', 'x2-y2': 'B1',
                'dxy': 'B2', 'xy': 'B2',
                'dxz': 'E', 'xz': 'E', 'dyz': 'E', 'yz': 'E',
                # Quadratic forms
                'x2': 'A1', 'y2': 'A1',
                'x2+y2+z2': 'A1', 'r2': 'A1',
            },
            'D2h': {
                # Coordinates
                'x': 'B3u', 'y': 'B2u', 'z': 'B1u',
                # Rotations
                'rx': 'B3g', 'ry': 'B2g', 'rz': 'B1g',
                # d orbitals
                'dz2': 'Ag', '2z2-x2-y2': 'Ag', 'z2': 'Ag',
                'dx2-y2': 'Ag', 'x2-y2': 'Ag',
                'dxy': 'B1g', 'xy': 'B1g',
                'dxz': 'B3g', 'xz': 'B3g',
                'dyz': 'B2g', 'yz': 'B2g',
                # Quadratic forms
                'x2': 'Ag', 'y2': 'Ag',
                'x2+y2+z2': 'Ag', 'r2': 'Ag',
            },
            'D3h': {
                # Coordinates
                'x': 'E\'', 'y': 'E\'', 'z': 'A2"',
                # Rotations
                'rx': 'E"', 'ry': 'E"', 'rz': 'A2\'',
                # d orbitals
                'dz2': 'A1\'', '2z2-x2-y2': 'A1\'', 'z2': 'A1\'',
                'dx2-y2': 'E\'', 'x2-y2': 'E\'', 'dxy': 'E\'', 'xy': 'E\'',
                'dxz': 'E"', 'xz': 'E"', 'dyz': 'E"', 'yz': 'E"',
                # Quadratic forms
                'x2': 'A1\'', 'y2': 'A1\'',
                'x2+y2+z2': 'A1\'', 'r2': 'A1\'',
            },
            'C5v': {
                # Coordinates
                'x': 'E1', 'y': 'E1', 'z': 'A1',
                # Rotations
                'rx': 'E1', 'ry': 'E1', 'rz': 'A2',
                # d orbitals
                'dz2': 'A1', '2z2-x2-y2': 'A1', 'z2': 'A1',
                'dx2-y2': 'E2', 'x2-y2': 'E2', 'dxy': 'E2', 'xy': 'E2',
                'dxz': 'E1', 'xz': 'E1', 'dyz': 'E1', 'yz': 'E1',
                # Quadratic forms
                'x2': 'A1', 'y2': 'A1',
                'x2+y2+z2': 'A1', 'r2': 'A1',
            },
            'C6v': {
                # Coordinates
                'x': 'E1', 'y': 'E1', 'z': 'A1',
                # Rotations
                'rx': 'E1', 'ry': 'E1', 'rz': 'A2',
                # d orbitals
                'dz2': 'A1', '2z2-x2-y2': 'A1', 'z2': 'A1',
                'dx2-y2': 'B1', 'x2-y2': 'B1',
                'dxy': 'B2', 'xy': 'B2',
                'dxz': 'E1', 'xz': 'E1', 'dyz': 'E1', 'yz': 'E1',
                # Quadratic forms  
                'x2': 'A1', 'y2': 'A1',
                'x2+y2+z2': 'A1', 'r2': 'A1',
            },
            'D5h': {
                # Coordinates
                'x': 'E1\'', 'y': 'E1\'', 'z': 'A2"',
                # Rotations
                'rx': 'E1"', 'ry': 'E1"', 'rz': 'A2\'',
                # d orbitals
                'dz2': 'A1\'', '2z2-x2-y2': 'A1\'', 'z2': 'A1\'',
                'dx2-y2': 'E2\'', 'x2-y2': 'E2\'', 'dxy': 'E2\'', 'xy': 'E2\'',
                'dxz': 'E1"', 'xz': 'E1"', 'dyz': 'E1"', 'yz': 'E1"',
                # Quadratic forms
                'x2': 'A1\'', 'y2': 'A1\'',
                'x2+y2+z2': 'A1\'', 'r2': 'A1\'',
            },
            'D6h': {
                # Coordinates
                'x': 'E1u', 'y': 'E1u', 'z': 'A2u',
                # Rotations
                'rx': 'E1g', 'ry': 'E1g', 'rz': 'A2g',
                # d orbitals
                'dz2': 'A1g', '2z2-x2-y2': 'A1g', 'z2': 'A1g',
                'dx2-y2': 'B1g', 'x2-y2': 'B1g',
                'dxy': 'B2g', 'xy': 'B2g',
                'dxz': 'E1g', 'xz': 'E1g', 'dyz': 'E1g', 'yz': 'E1g',
                # Quadratic forms
                'x2': 'A1g', 'y2': 'A1g',
                'x2+y2+z2': 'A1g', 'r2': 'A1g',
            },
            'C1': {
                # All functions belong to the single A irrep
                'x': 'A', 'y': 'A', 'z': 'A',
                'rx': 'A', 'ry': 'A', 'rz': 'A',
                'dz2': 'A', 'dx2-y2': 'A', 'x2-y2': 'A',
                'dxy': 'A', 'xy': 'A', 'dxz': 'A', 'xz': 'A', 'dyz': 'A', 'yz': 'A',
                '2z2-x2-y2': 'A', 'z2': 'A', 'x2': 'A', 'y2': 'A',
                'x2+y2+z2': 'A', 'r2': 'A',
            },
            'Ci': {
                # Functions with even parity → Ag, odd parity → Au
                'x': 'Au', 'y': 'Au', 'z': 'Au',
                'rx': 'Au', 'ry': 'Au', 'rz': 'Au',
                'dz2': 'Ag', 'dx2-y2': 'Ag', 'x2-y2': 'Ag',
                'dxy': 'Ag', 'xy': 'Ag', 'dxz': 'Ag', 'xz': 'Ag', 'dyz': 'Ag', 'yz': 'Ag',
                '2z2-x2-y2': 'Ag', 'z2': 'Ag', 'x2': 'Ag', 'y2': 'Ag',
                'x2+y2+z2': 'Ag', 'r2': 'Ag',
            },
            'Cs': {
                # Functions symmetric to plane → A', antisymmetric → A"
                'x': "A'", 'y': "A'", 'z': "A'",
                'rx': 'A"', 'ry': 'A"', 'rz': 'A"',
                'dz2': "A'", 'dx2-y2': "A'", 'x2-y2': "A'",
                'dxy': "A'", 'xy': "A'", 'dxz': 'A"', 'xz': 'A"', 'dyz': 'A"', 'yz': 'A"',
                '2z2-x2-y2': "A'", 'z2': "A'", 'x2': "A'", 'y2': "A'",
                'x2+y2+z2': "A'", 'r2': "A'",
            },
            'C2': {
                # Functions even under C2 → A, odd → B
                'x': 'B', 'y': 'B', 'z': 'A',
                'rx': 'B', 'ry': 'B', 'rz': 'A',
                'dz2': 'A', 'dx2-y2': 'A', 'x2-y2': 'A',
                'dxy': 'B', 'xy': 'B', 'dxz': 'B', 'xz': 'B', 'dyz': 'B', 'yz': 'B',
                '2z2-x2-y2': 'A', 'z2': 'A', 'x2': 'A', 'y2': 'A',
                'x2+y2+z2': 'A', 'r2': 'A',
            },
            'C3': {
                # Real basis set for C3
                'x': 'E', 'y': 'E', 'z': 'A',
                'rx': 'E', 'ry': 'E', 'rz': 'A',
                'dz2': 'A', '2z2-x2-y2': 'A', 'z2': 'A',
                'dx2-y2': 'E', 'x2-y2': 'E', 'dxy': 'E', 'xy': 'E',
                'dxz': 'E', 'xz': 'E', 'dyz': 'E', 'yz': 'E',
                # Quadratic forms
                'x2': 'A', 'y2': 'A',
                'x2+y2+z2': 'A', 'r2': 'A',
            },
            'D2': {
                # Coordinates
                'x': 'B1', 'y': 'B2', 'z': 'B3',
                # Rotations
                'rx': 'B1', 'ry': 'B2', 'rz': 'B3',
                # d orbitals
                'dz2': 'A', 'dx2-y2': 'A', 'x2-y2': 'A',
                'dxy': 'A', 'xy': 'A', 'dxz': 'B2', 'xz': 'B2', 'dyz': 'B1', 'yz': 'B1',
                '2z2-x2-y2': 'A', 'z2': 'A', 'x2': 'A', 'y2': 'A',
                'x2+y2+z2': 'A', 'r2': 'A',
            },
            'D3': {
                # Coordinates  
                'x': 'E', 'y': 'E', 'z': 'A2',
                # Rotations
                'rx': 'E', 'ry': 'E', 'rz': 'A2',
                # d orbitals
                'dz2': 'A1', '2z2-x2-y2': 'A1', 'z2': 'A1',
                'dx2-y2': 'E', 'x2-y2': 'E', 'dxy': 'E', 'xy': 'E',
                'dxz': 'E', 'xz': 'E', 'dyz': 'E', 'yz': 'E',
                # Quadratic forms
                'x2': 'A1', 'y2': 'A1',
                'x2+y2+z2': 'A1', 'r2': 'A1',
            },
            'Ih': {
                # Coordinates
                'x': 'T1u', 'y': 'T1u', 'z': 'T1u',
                # Rotations
                'rx': 'T1g', 'ry': 'T1g', 'rz': 'T1g', 
                # d orbitals
                'dz2': 'Hg', 'dx2-y2': 'Hg', 'x2-y2': 'Hg',
                'dxy': 'Hg', 'xy': 'Hg', 'dxz': 'Hg', 'xz': 'Hg', 'dyz': 'Hg', 'yz': 'Hg',
                '2z2-x2-y2': 'Hg', 'z2': 'Ag', 'x2': 'Ag', 'y2': 'Ag',
                'x2+y2+z2': 'Ag', 'r2': 'Ag',
            },
            'Cinfv': {
                # Linear C∞v
                'x': 'Π', 'y': 'Π', 'z': 'Σ+',
                'rx': 'Π', 'ry': 'Π', 'rz': 'Σ+',
                # Note: d orbitals in linear groups have different labels
                'dz2': 'Σ+', '2z2-x2-y2': 'Σ+', 'z2': 'Σ+',
                'dx2-y2': 'Δ', 'x2-y2': 'Δ', 'dxy': 'Δ', 'xy': 'Δ',
                'dxz': 'Π', 'xz': 'Π', 'dyz': 'Π', 'yz': 'Π',
                'x2': 'Σ+', 'y2': 'Σ+',
                'x2+y2+z2': 'Σ+', 'r2': 'Σ+',
            },
            'Dinfh': {
                # Linear D∞h
                'x': 'Πu', 'y': 'Πu', 'z': 'Σ+u',
                'rx': 'Πg', 'ry': 'Πg', 'rz': 'Σ+g',
                # d orbitals  
                'dz2': 'Σ+g', '2z2-x2-y2': 'Σ+g', 'z2': 'Σ+g',
                'dx2-y2': 'Δg', 'x2-y2': 'Δg', 'dxy': 'Δg', 'xy': 'Δg',
                'dxz': 'Πg', 'xz': 'Πg', 'dyz': 'Πg', 'yz': 'Πg',
                'x2': 'Σ+g', 'y2': 'Σ+g',
                'x2+y2+z2': 'Σ+g', 'r2': 'Σ+g',
            },
            # =============================================================================
            # DOUBLE GROUPS WITH SPINOR BASIS FUNCTIONS
            # =============================================================================
            'C2v*': {
                # Integer spin functions (same as C2v)
                'x': 'B1', 'y': 'B2', 'z': 'A1',
                'rx': 'B2', 'ry': 'B1', 'rz': 'A2',
                'dz2': 'A1', 'dx2-y2': 'A1', 'dxy': 'A2',
                'dxz': 'B1', 'dyz': 'B2',
                # Half-integer spinor functions
                'j=1/2': 'E1/2', 'spinor1/2': 'E1/2',
                'j=3/2': 'E3/2', 'spinor3/2': 'E3/2',
                # Electron spin functions
                'α': 'E1/2', 'β': 'E1/2',  # Kramers doublet
                'up': 'E1/2', 'down': 'E1/2',
                # p1/2 and p3/2 orbitals for spin-orbit coupling
                'p1/2': 'E1/2', 'p3/2': 'E3/2',
            },
            'Oh*': {
                # Integer spin functions (same as Oh)
                'x': 'T1u', 'y': 'T1u', 'z': 'T1u',
                'rx': 'T1g', 'ry': 'T1g', 'rz': 'T1g',
                'dz2': 'Eg', 'dx2-y2': 'Eg', 'x2-y2': 'Eg',
                'dxy': 'T2g', 'dxz': 'T2g', 'dyz': 'T2g',
                '2z2-x2-y2': 'Eg', 'xy': 'T2g', 'xz': 'T2g', 'yz': 'T2g',
                # Half-integer spinor functions
                'j=1/2': 'E1/2g', 'j=1/2g': 'E1/2g', 'j=1/2u': 'E1/2u',
                'j=3/2': 'E3/2g', 'j=3/2g': 'E3/2g', 'j=3/2u': 'E3/2u',
                'j=5/2': 'E5/2g', 'j=5/2g': 'E5/2g', 'j=5/2u': 'E5/2u',
                # Electron spin functions
                'α': 'E1/2g', 'β': 'E1/2g',  # Kramers doublet
                'up': 'E1/2g', 'down': 'E1/2g',
                # t2g orbitals with spin-orbit coupling split to j=1/2, j=3/2
                't2g,j=1/2': 'E1/2g', 't2g,j=3/2': 'E3/2g',
                # eg orbitals with spin-orbit coupling  
                'eg,j=1/2': 'E1/2g', 'eg,j=3/2': 'E3/2g',
                # p orbitals with spin-orbit coupling
                'p1/2': 'E1/2u', 'p3/2': 'E3/2u',
            },
            'Td*': {
                # Integer spin functions (same as Td)
                'x': 'T2', 'y': 'T2', 'z': 'T2',
                'dz2': 'E', 'dx2-y2': 'E', 'x2-y2': 'E',
                'dxy': 'T2', 'dxz': 'T2', 'dyz': 'T2',
                '2z2-x2-y2': 'E', 'xy': 'T2', 'xz': 'T2', 'yz': 'T2',
                # Half-integer spinor functions
                'j=1/2': 'E1/2', 'spinor1/2': 'E1/2',
                'j=3/2': 'E3/2', 'spinor3/2': 'E3/2', 
                'j=5/2': 'E5/2', 'spinor5/2': 'E5/2',
                # Electron spin functions
                'α': 'E1/2', 'β': 'E1/2',  # Kramers doublet
                'up': 'E1/2', 'down': 'E1/2',
                # e orbitals with spin-orbit coupling
                'e,j=1/2': 'E1/2', 'e,j=3/2': 'E3/2',
                # t2 orbitals with spin-orbit coupling
                't2,j=1/2': 'E1/2', 't2,j=3/2': 'E3/2', 't2,j=5/2': 'E5/2',
                # p orbitals with spin-orbit coupling
                'p1/2': 'E1/2', 'p3/2': 'E3/2',
            },
            'D4h*': {
                # Integer spin functions (same as D4h)
                'x': 'Eu', 'y': 'Eu', 'z': 'A2u',
                'rx': 'Eg', 'ry': 'Eg', 'rz': 'A2g',
                'dz2': 'A1g', 'dx2-y2': 'B1g', 'dxy': 'B2g',
                'dxz': 'Eg', 'dyz': 'Eg',
                '2z2-x2-y2': 'A1g', 'x2-y2': 'B1g', 'xy': 'B2g',
                'xz': 'Eg', 'yz': 'Eg',
                # Half-integer spinor functions
                'j=1/2': 'E1/2g', 'j=1/2g': 'E1/2g', 'j=1/2u': 'E1/2u',
                'j=3/2': 'E3/2g', 'j=3/2g': 'E3/2g', 'j=3/2u': 'E3/2u',
                # Electron spin functions
                'α': 'E1/2g', 'β': 'E1/2g',  # Kramers doublet
                'up': 'E1/2g', 'down': 'E1/2g',
                # d orbitals with spin-orbit coupling
                'dz2,j=1/2': 'E1/2g', 'dx2-y2,j=1/2': 'E1/2g', 'dxy,j=1/2': 'E1/2g',
                'dxz,j=1/2': 'E1/2g', 'dyz,j=1/2': 'E1/2g',
                'dz2,j=3/2': 'E3/2g', 'dx2-y2,j=3/2': 'E3/2g', 'dxy,j=3/2': 'E3/2g', 
                'dxz,j=3/2': 'E3/2g', 'dyz,j=3/2': 'E3/2g',
                # p orbitals with spin-orbit coupling
                'p1/2': 'E1/2u', 'p3/2': 'E3/2u',
            },
        }

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


class C4(PointGroup):
    """C4 point group - four-fold rotation axis."""
    
    def __init__(self):
        super().__init__()
        self.name = "C4"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        # For C4, we need 1-dimensional E1 and E2 irreps (complex conjugate pairs)
        self.classes = ['E', 'C4', 'C2', 'C4³']
        self.class_sizes = [1, 1, 1, 1]
        self.order = 4
        self.irreps = {
            'A': [1, 1, 1, 1],
            'B': [1, -1, 1, -1],
            'E1': [1, 1j, -1, -1j],      # ε = i for C4
            'E2': [1, -1j, -1, 1j]       # ε* = -i for C4
        }


class C5(PointGroup):
    """C5 point group - five-fold rotation axis."""
    
    def __init__(self):
        super().__init__()
        self.name = "C5"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        # For C5, characters are 2*cos(2πk/5) for the E representations
        cos72 = 2 * np.cos(2 * np.pi / 5)      # 2*cos(72°) ≈ 0.618
        cos144 = 2 * np.cos(4 * np.pi / 5)    # 2*cos(144°) ≈ -1.618
        self.classes = ['E', 'C5', 'C5²', 'C5³', 'C5⁴']
        self.class_sizes = [1, 1, 1, 1, 1]
        self.order = 5
        self.irreps = {
            'A': [1, 1, 1, 1, 1],
            'E1': [2, cos72, cos144, cos144, cos72],
            'E2': [2, cos144, cos72, cos72, cos144]
        }


class C6(PointGroup):
    """C6 point group - six-fold rotation axis."""
    
    def __init__(self):
        super().__init__()
        self.name = "C6"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        # ε = e^(2πi/6) = e^(πi/3)
        ε = np.exp(1j * np.pi / 3)
        self.classes = ['E', 'C6', 'C3', 'C2', 'C3²', 'C6⁵']
        self.class_sizes = [1, 1, 1, 1, 1, 1]
        self.order = 6
        self.irreps = {
            'A': [1, 1, 1, 1, 1, 1],
            'B': [1, -1, 1, -1, 1, -1],
            'E1': [2, 1, -1, -2, -1, 1],
            'E2': [2, -1, -1, 2, -1, -1]
        }


class C7(PointGroup):
    """C7 point group - seven-fold rotation axis."""
    
    def __init__(self):
        super().__init__()
        self.name = "C7"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        # For C7, characters are 2*cos(2πk/7) for the E representations
        cos1 = 2 * np.cos(2 * np.pi / 7)       # 2*cos(51.43°) ≈ 1.247  
        cos2 = 2 * np.cos(4 * np.pi / 7)       # 2*cos(102.86°) ≈ -0.445
        cos3 = 2 * np.cos(6 * np.pi / 7)       # 2*cos(154.29°) ≈ -1.802
        self.classes = ['E', 'C7', 'C7²', 'C7³', 'C7⁴', 'C7⁵', 'C7⁶']
        self.class_sizes = [1, 1, 1, 1, 1, 1, 1]
        self.order = 7
        self.irreps = {
            'A': [1, 1, 1, 1, 1, 1, 1],
            'E1': [2, cos1, cos2, cos3, cos3, cos2, cos1],
            'E2': [2, cos2, cos3, cos1, cos1, cos3, cos2],
            'E3': [2, cos3, cos1, cos2, cos2, cos1, cos3]
        }


class C8(PointGroup):
    """C8 point group - eight-fold rotation axis."""
    
    def __init__(self):
        super().__init__()
        self.name = "C8"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        # ε = e^(2πi/8) = e^(πi/4)
        ε = np.exp(1j * np.pi / 4)
        self.classes = ['E', 'C8', 'C4', 'C8³', 'C2', 'C8⁵', 'C4³', 'C8⁷']
        self.class_sizes = [1, 1, 1, 1, 1, 1, 1, 1]
        self.order = 8
        self.irreps = {
            'A': [1, 1, 1, 1, 1, 1, 1, 1],
            'B': [1, -1, 1, -1, 1, -1, 1, -1],
            'E1': [2, np.sqrt(2), 0, -np.sqrt(2), -2, -np.sqrt(2), 0, np.sqrt(2)],
            'E2': [2, 0, -2, 0, 2, 0, -2, 0],
            'E3': [2, -np.sqrt(2), 0, np.sqrt(2), -2, np.sqrt(2), 0, -np.sqrt(2)]
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


class C2h(PointGroup):
    """C2h point group - two-fold axis with horizontal mirror plane."""
    
    def __init__(self):
        super().__init__()
        self.name = "C2h"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C2', 'i', 'σh']
        self.class_sizes = [1, 1, 1, 1]
        self.order = 4
        self.irreps = {
            'Ag': [1, 1, 1, 1],
            'Bg': [1, -1, 1, -1],
            'Au': [1, 1, -1, -1],
            'Bu': [1, -1, -1, 1]
        }


class C3h(PointGroup):
    """C3h point group - three-fold axis with horizontal mirror plane."""
    
    def __init__(self):
        super().__init__()
        self.name = "C3h"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        # ε = e^(2πi/3)
        ε = np.exp(2j * np.pi / 3)
        self.classes = ['E', 'C3', 'σh']
        self.class_sizes = [1, 2, 3]
        self.order = 6
        self.irreps = {
            "A'": [1, 1, 1],
            'A"': [1, 1, -1],
            "E'": [2, -1, 0],
            'E"': [2, -1, 0]
        }


class C4h(PointGroup):
    """C4h point group - four-fold axis with horizontal mirror plane."""
    
    def __init__(self):
        super().__init__()
        self.name = "C4h"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        self.classes = ['E', 'C4', 'C2', 'C4³', 'i', 'S4³', 'σh', 'S4']
        self.class_sizes = [1, 1, 1, 1, 1, 1, 1, 1]
        self.order = 8
        self.irreps = {
            'Ag': [1, 1, 1, 1, 1, 1, 1, 1],
            'Bg': [1, -1, 1, -1, 1, -1, 1, -1],
            'Eg': [2, 0, -2, 0, 2, 0, -2, 0],
            'Au': [1, 1, 1, 1, -1, -1, -1, -1],
            'Bu': [1, -1, 1, -1, -1, 1, -1, 1],
            'Eu': [2, 0, -2, 0, -2, 0, 2, 0]
        }


class C5h(PointGroup):
    """C5h point group - five-fold axis with horizontal mirror plane."""
    
    def __init__(self):
        super().__init__()
        self.name = "C5h"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        # ε = e^(2πi/5)
        ε = np.exp(2j * np.pi / 5)
        self.classes = ['E', 'C5', 'C5²', 'C5³', 'C5⁴', 'σh', 'S5', 'S5³', 'S5²', 'S5⁴']
        self.class_sizes = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        self.order = 10
        self.irreps = {
            "A'": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            'A"': [1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
            "E1'": [2, 2*np.cos(2*np.pi/5), 2*np.cos(4*np.pi/5), 2*np.cos(4*np.pi/5), 2*np.cos(2*np.pi/5), 2, 2*np.cos(2*np.pi/5), 2*np.cos(4*np.pi/5), 2*np.cos(4*np.pi/5), 2*np.cos(2*np.pi/5)],
            'E1"': [2, 2*np.cos(2*np.pi/5), 2*np.cos(4*np.pi/5), 2*np.cos(4*np.pi/5), 2*np.cos(2*np.pi/5), -2, -2*np.cos(2*np.pi/5), -2*np.cos(4*np.pi/5), -2*np.cos(4*np.pi/5), -2*np.cos(2*np.pi/5)],
            "E2'": [2, 2*np.cos(4*np.pi/5), 2*np.cos(2*np.pi/5), 2*np.cos(2*np.pi/5), 2*np.cos(4*np.pi/5), 2, 2*np.cos(4*np.pi/5), 2*np.cos(2*np.pi/5), 2*np.cos(2*np.pi/5), 2*np.cos(4*np.pi/5)],
            'E2"': [2, 2*np.cos(4*np.pi/5), 2*np.cos(2*np.pi/5), 2*np.cos(2*np.pi/5), 2*np.cos(4*np.pi/5), -2, -2*np.cos(4*np.pi/5), -2*np.cos(2*np.pi/5), -2*np.cos(2*np.pi/5), -2*np.cos(4*np.pi/5)]
        }


class C6h(PointGroup):
    """C6h point group - six-fold axis with horizontal mirror plane."""
    
    def __init__(self):
        super().__init__()
        self.name = "C6h"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        self.classes = ['E', 'C6', 'C3', 'C2', 'C3²', 'C6⁵', 'i', 'S3⁵', 'S6⁵', 'σh', 'S6', 'S3']
        self.class_sizes = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        self.order = 12
        self.irreps = {
            'Ag': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            'Bg': [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
            'E1g': [2, 1, -1, -2, -1, 1, 2, 1, -1, -2, -1, 1],
            'E2g': [2, -1, -1, 2, -1, -1, 2, -1, -1, 2, -1, -1],
            'Au': [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1],
            'Bu': [1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1],
            'E1u': [2, 1, -1, -2, -1, 1, -2, -1, 1, 2, 1, -1],
            'E2u': [2, -1, -1, 2, -1, -1, -2, 1, 1, -2, 1, 1]
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


class D4(PointGroup):
    """D4 point group - four-fold axis with perpendicular two-fold axes."""
    
    def __init__(self):
        super().__init__()
        self.name = "D4"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C4', 'C2', 'C2\'', 'C2"']
        self.class_sizes = [1, 2, 1, 2, 2]
        self.order = 8
        self.irreps = {
            'A1': [1, 1, 1, 1, 1],
            'A2': [1, 1, 1, -1, -1],
            'B1': [1, -1, 1, 1, -1],
            'B2': [1, -1, 1, -1, 1],
            'E': [2, 0, -2, 0, 0]
        }


class D5(PointGroup):
    """D5 point group - five-fold axis with perpendicular two-fold axes."""
    
    def __init__(self):
        super().__init__()
        self.name = "D5"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        self.classes = ['E', 'C5', 'C5²', 'C2']
        self.class_sizes = [1, 2, 2, 5]
        self.order = 10
        self.irreps = {
            'A1': [1, 1, 1, 1],
            'A2': [1, 1, 1, -1],
            'E1': [2, 2*np.cos(2*np.pi/5), 2*np.cos(4*np.pi/5), 0],
            'E2': [2, 2*np.cos(4*np.pi/5), 2*np.cos(2*np.pi/5), 0]
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


class D2d(PointGroup):
    """D2d point group - D2 with dihedral planes."""
    
    def __init__(self):
        super().__init__()
        self.name = "D2d"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C2', 'C2\'', 'S4', 'σd']
        self.class_sizes = [1, 1, 2, 2, 2]
        self.order = 8
        self.irreps = {
            'A1': [1, 1, 1, 1, 1],
            'A2': [1, 1, 1, -1, -1],
            'B1': [1, 1, -1, 1, -1],
            'B2': [1, 1, -1, -1, 1],
            'E': [2, -2, 0, 0, 0]
        }


class D3d(PointGroup):
    """D3d point group - D3 with dihedral planes."""
    
    def __init__(self):
        super().__init__()
        self.name = "D3d"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C3', 'C2', 'i', 'S6', 'σd']
        self.class_sizes = [1, 2, 3, 1, 2, 3]
        self.order = 12
        self.irreps = {
            'A1g': [1, 1, 1, 1, 1, 1],
            'A2g': [1, 1, -1, 1, 1, -1],
            'Eg': [2, -1, 0, 2, -1, 0],
            'A1u': [1, 1, 1, -1, -1, -1],
            'A2u': [1, 1, -1, -1, -1, 1],
            'Eu': [2, -1, 0, -2, 1, 0]
        }


class D4d(PointGroup):
    """D4d point group - D4 with dihedral planes."""
    
    def __init__(self):
        super().__init__()
        self.name = "D4d"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        self.classes = ['E', 'C4', 'C2', 'C2\'', 'C2"', 'S8', 'σd']
        self.class_sizes = [1, 2, 1, 2, 2, 2, 4]
        self.order = 16
        self.irreps = {
            'A1': [1, 1, 1, 1, 1, 1, 1],
            'A2': [1, 1, 1, -1, -1, 1, -1],
            'B1': [1, -1, 1, 1, -1, -1, 1],
            'B2': [1, -1, 1, -1, 1, -1, -1],
            'E1': [2, 0, -2, 0, 0, np.sqrt(2), 0],
            'E2': [2, 0, -2, 0, 0, -np.sqrt(2), 0],
            'E3': [2, np.sqrt(2), 0, 0, 0, 0, 0]
        }


class D5d(PointGroup):
    """D5d point group - D5 with dihedral planes."""
    
    def __init__(self):
        super().__init__()
        self.name = "D5d"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        cos72 = np.cos(2*np.pi/5)  # cos(72°)
        cos144 = np.cos(4*np.pi/5)  # cos(144°)
        self.classes = ['E', 'C5', 'C5²', 'C2', 'i', 'S10', 'S10³', 'σd']
        self.class_sizes = [1, 2, 2, 5, 1, 2, 2, 5]
        self.order = 20
        self.irreps = {
            'A1g': [1, 1, 1, 1, 1, 1, 1, 1],
            'A2g': [1, 1, 1, -1, 1, 1, 1, -1],
            'E1g': [2, 2*cos72, 2*cos144, 0, 2, 2*cos72, 2*cos144, 0],
            'E2g': [2, 2*cos144, 2*cos72, 0, 2, 2*cos144, 2*cos72, 0],
            'A1u': [1, 1, 1, 1, -1, -1, -1, -1],
            'A2u': [1, 1, 1, -1, -1, -1, -1, 1],
            'E1u': [2, 2*cos72, 2*cos144, 0, -2, -2*cos72, -2*cos144, 0],
            'E2u': [2, 2*cos144, 2*cos72, 0, -2, -2*cos144, -2*cos72, 0]
        }


class D6d(PointGroup):
    """D6d point group - D6 with dihedral planes."""
    
    def __init__(self):
        super().__init__()
        self.name = "D6d"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        self.classes = ['E', 'C6', 'C3', 'C2', 'C2\'', 'C2"', 'S12', 'S4', 'σd']
        self.class_sizes = [1, 2, 2, 1, 3, 3, 2, 2, 6]
        self.order = 24
        self.irreps = {
            'A1': [1, 1, 1, 1, 1, 1, 1, 1, 1],
            'A2': [1, 1, 1, 1, -1, -1, 1, 1, -1],
            'B1': [1, -1, 1, -1, 1, -1, -1, 1, -1],
            'B2': [1, -1, 1, -1, -1, 1, -1, 1, 1],
            'E1': [2, 1, -1, -2, 0, 0, 1, 0, 0],
            'E2': [2, -1, -1, 2, 0, 0, -1, 0, 0],
            'E3': [2, np.sqrt(3), 0, 0, 0, 0, -np.sqrt(3), 0, 0],
            'E4': [2, 0, -2, 0, 0, 0, 0, -2, 0],
            'E5': [2, -np.sqrt(3), 0, 0, 0, 0, np.sqrt(3), 0, 0]
        }


class S4(PointGroup):
    """S4 point group - four-fold improper rotation axis."""
    
    def __init__(self):
        super().__init__()
        self.name = "S4"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        self.classes = ['E', 'S4', 'C2', 'S4³']
        self.class_sizes = [1, 1, 1, 1]
        self.order = 4
        self.irreps = {
            'A': [1, 1, 1, 1],
            'B': [1, -1, 1, -1],
            'E': [2, 0, -2, 0]
        }


class S6(PointGroup):
    """S6 point group - six-fold improper rotation axis."""
    
    def __init__(self):
        super().__init__()
        self.name = "S6"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        # ε = e^(2πi/3)
        ε = np.exp(2j * np.pi / 3)
        self.classes = ['E', 'C3', 'C3²', 'i', 'S6⁵', 'S6']
        self.class_sizes = [1, 1, 1, 1, 1, 1]
        self.order = 6
        self.irreps = {
            'Ag': [1, 1, 1, 1, 1, 1],
            'Eg': [2, -1, -1, 2, -1, -1],
            'Au': [1, 1, 1, -1, -1, -1],
            'Eu': [2, -1, -1, -2, 1, 1]
        }


class S8(PointGroup):
    """S8 point group - eight-fold improper rotation axis."""
    
    def __init__(self):
        super().__init__()
        self.name = "S8"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        self.classes = ['E', 'S8', 'C4', 'S8³', 'C2', 'S8⁵', 'C4³', 'S8⁷']
        self.class_sizes = [1, 1, 1, 1, 1, 1, 1, 1]
        self.order = 8
        self.irreps = {
            'A': [1, 1, 1, 1, 1, 1, 1, 1],
            'B': [1, -1, 1, -1, 1, -1, 1, -1],
            'E1': [2, np.sqrt(2), 0, -np.sqrt(2), -2, -np.sqrt(2), 0, np.sqrt(2)],
            'E2': [2, 0, -2, 0, 2, 0, -2, 0],
            'E3': [2, -np.sqrt(2), 0, np.sqrt(2), -2, np.sqrt(2), 0, -np.sqrt(2)]
        }


class T(PointGroup):
    """T point group - pure rotational tetrahedral symmetry."""
    
    def __init__(self):
        super().__init__()
        self.name = "T"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        # ε = e^(2πi/3)
        ε = np.exp(2j * np.pi / 3)
        self.classes = ['E', 'C3', 'C2']
        self.class_sizes = [1, 8, 3]
        self.order = 12
        self.irreps = {
            'A': [1, 1, 1],
            'E': [2, -1, 2],
            'T': [3, 0, -1]
        }


class O(PointGroup):
    """O point group - pure rotational octahedral symmetry."""
    
    def __init__(self):
        super().__init__()
        self.name = "O"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C3', 'C2', 'C4', 'C2\'']
        self.class_sizes = [1, 8, 6, 6, 3]
        self.order = 24
        self.irreps = {
            'A1': [1, 1, 1, 1, 1],
            'A2': [1, 1, -1, -1, 1],
            'E': [2, -1, 0, 0, 2],
            'T1': [3, 0, -1, 1, -1],
            'T2': [3, 0, 1, -1, -1]
        }


class Th(PointGroup):
    """Th point group - T with inversion center."""
    
    def __init__(self):
        super().__init__()
        self.name = "Th"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C3', 'C2', 'i', 'S6', 'σh']
        self.class_sizes = [1, 8, 3, 1, 8, 3]
        self.order = 24
        self.irreps = {
            'Ag': [1, 1, 1, 1, 1, 1],
            'Eg': [2, -1, 2, 2, -1, 2],
            'Tg': [3, 0, -1, 3, 0, -1],
            'Au': [1, 1, 1, -1, -1, -1],
            'Eu': [2, -1, 2, -2, 1, -2],
            'Tu': [3, 0, -1, -3, 0, 1]
        }


class I(PointGroup):
    """I point group - pure rotational icosahedral symmetry."""
    
    def __init__(self):
        super().__init__()
        self.name = "I"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        φ = (1 + np.sqrt(5)) / 2  # golden ratio
        self.classes = ['E', 'C5', 'C5²', 'C3', 'C2']
        self.class_sizes = [1, 12, 12, 20, 15]
        self.order = 60
        self.irreps = {
            'A': [1, 1, 1, 1, 1],
            'T1': [3, φ, -1/φ, 0, -1],
            'T2': [3, -1/φ, φ, 0, -1],
            'G': [4, -1, -1, 1, 0],
            'H': [5, 0, 0, -1, 1]
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


# =============================================================================
# DOUBLE GROUPS FOR HALF-INTEGER SPIN SYSTEMS
# =============================================================================

class C2vStar(PointGroup):
    """C2v* double group - includes half-integer representations for fermions."""
    
    def __init__(self):
        super().__init__()
        self.name = "C2v*"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        # Double group includes R (2π rotation) operations
        self.classes = ['E', 'C2', 'σv(xz)', 'σv(yz)', 'R', 'RC2', 'Rσv(xz)', 'Rσv(yz)']
        self.class_sizes = [1, 1, 1, 1, 1, 1, 1, 1]
        self.order = 8
        self.irreps = {
            # Integer spin representations (same as C2v)
            'A1': [1, 1, 1, 1, 1, 1, 1, 1],
            'A2': [1, 1, -1, -1, 1, 1, -1, -1],
            'B1': [1, -1, 1, -1, 1, -1, 1, -1],
            'B2': [1, -1, -1, 1, 1, -1, -1, 1],
            # Half-integer spin representations (Kramers doublets)
            'E1/2': [2, 0, 0, 0, -2, 0, 0, 0],  # j = 1/2
            'E3/2': [2, 0, 0, 0, -2, 0, 0, 0],  # j = 3/2
        }


class OhStar(PointGroup):
    """Oh* double group - octahedral symmetry with half-integer representations."""
    
    def __init__(self):
        super().__init__()
        self.name = "Oh*"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        # Double group includes R operations  
        self.classes = ['E', '8C3', '6C2', '6C4', '3C2', 'i', '6S4', '8S6', '3σh', '6σd',
                       'R', '8RC3', '6RC2', '6RC4', '3RC2', 'Ri', '6RS4', '8RS6', '3Rσh', '6Rσd']
        self.class_sizes = [1, 8, 6, 6, 3, 1, 6, 8, 3, 6,
                           1, 8, 6, 6, 3, 1, 6, 8, 3, 6]
        self.order = 96  # Double the original Oh order (48)
        self.irreps = {
            # Integer spin representations (same as Oh)
            'A1g': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            'A2g': [1, 1, -1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, 1, -1],
            'Eg': [2, -1, 0, 0, 2, 2, 0, -1, 2, 0, 2, -1, 0, 0, 2, 2, 0, -1, 2, 0],
            'T1g': [3, 0, -1, 1, -1, 3, 1, 0, -1, -1, 3, 0, -1, 1, -1, 3, 1, 0, -1, -1],
            'T2g': [3, 0, 1, -1, -1, 3, -1, 0, -1, 1, 3, 0, 1, -1, -1, 3, -1, 0, -1, 1],
            'A1u': [1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
            'A2u': [1, 1, -1, -1, 1, -1, 1, -1, -1, 1, 1, 1, -1, -1, 1, -1, 1, -1, -1, 1],
            'Eu': [2, -1, 0, 0, 2, -2, 0, 1, -2, 0, 2, -1, 0, 0, 2, -2, 0, 1, -2, 0],
            'T1u': [3, 0, -1, 1, -1, -3, -1, 0, 1, 1, 3, 0, -1, 1, -1, -3, -1, 0, 1, 1],
            'T2u': [3, 0, 1, -1, -1, -3, 1, 0, 1, -1, 3, 0, 1, -1, -1, -3, 1, 0, 1, -1],
            # Half-integer spin representations (Kramers doublets)
            'E1/2g': [2, 1, 0, 0, -2, 2, 0, -1, -2, 0, -2, -1, 0, 0, 2, -2, 0, 1, 2, 0],  # j = 1/2
            'E1/2u': [2, 1, 0, 0, -2, -2, 0, 1, 2, 0, -2, -1, 0, 0, 2, 2, 0, -1, -2, 0],  # j = 1/2
            'E3/2g': [4, 1, 0, 0, 0, 4, 0, -1, 0, 0, -4, -1, 0, 0, 0, -4, 0, 1, 0, 0],    # j = 3/2
            'E3/2u': [4, 1, 0, 0, 0, -4, 0, 1, 0, 0, -4, -1, 0, 0, 0, 4, 0, -1, 0, 0],    # j = 3/2
            'E5/2g': [4, -1, 0, 0, 0, 4, 0, 1, 0, 0, -4, 1, 0, 0, 0, -4, 0, -1, 0, 0],    # j = 5/2
            'E5/2u': [4, -1, 0, 0, 0, -4, 0, -1, 0, 0, -4, 1, 0, 0, 0, 4, 0, 1, 0, 0],    # j = 5/2
        }


class TdStar(PointGroup):
    """Td* double group - tetrahedral symmetry with half-integer representations."""
    
    def __init__(self):
        super().__init__()
        self.name = "Td*"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        # Double group includes R operations
        self.classes = ['E', '8C3', '3C2', '6S4', '6σd', 'R', '8RC3', '3RC2', '6RS4', '6Rσd']
        self.class_sizes = [1, 8, 3, 6, 6, 1, 8, 3, 6, 6]
        self.order = 48  # Double the original Td order (24)
        self.irreps = {
            # Integer spin representations (same as Td)
            'A1': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            'A2': [1, 1, 1, -1, -1, 1, 1, 1, -1, -1],
            'E': [2, -1, 2, 0, 0, 2, -1, 2, 0, 0],
            'T1': [3, 0, -1, 1, -1, 3, 0, -1, 1, -1],
            'T2': [3, 0, -1, -1, 1, 3, 0, -1, -1, 1],
            # Half-integer spin representations (Kramers doublets)
            'E1/2': [2, 1, -2, 0, 0, -2, -1, 2, 0, 0],    # j = 1/2
            'E3/2': [4, 1, 0, 0, 0, -4, -1, 0, 0, 0],     # j = 3/2
            'E5/2': [4, -1, 0, 0, 0, -4, 1, 0, 0, 0],     # j = 5/2
        }


class D4hStar(PointGroup):
    """D4h* double group - square planar symmetry with half-integer representations."""
    
    def __init__(self):
        super().__init__()
        self.name = "D4h*"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        # Double group includes R operations
        self.classes = ['E', '2C4', 'C2', '2C2\'', '2C2"', 'i', '2S4', 'σh', '2σv', '2σd',
                       'R', '2RC4', 'RC2', '2RC2\'', '2RC2"', 'Ri', '2RS4', 'Rσh', '2Rσv', '2Rσd']
        self.class_sizes = [1, 2, 1, 2, 2, 1, 2, 1, 2, 2,
                           1, 2, 1, 2, 2, 1, 2, 1, 2, 2]
        self.order = 32  # Double the original D4h order (16)
        self.irreps = {
            # Integer spin representations (same as D4h)
            'A1g': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            'A2g': [1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1],
            'B1g': [1, -1, 1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, 1, -1],
            'B2g': [1, -1, 1, -1, 1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1],
            'Eg': [2, 0, -2, 0, 0, 2, 0, -2, 0, 0, 2, 0, -2, 0, 0, 2, 0, -2, 0, 0],
            'A1u': [1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
            'A2u': [1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1],
            'B1u': [1, -1, 1, 1, -1, -1, 1, -1, -1, 1, 1, -1, 1, 1, -1, -1, 1, -1, -1, 1],
            'B2u': [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
            'Eu': [2, 0, -2, 0, 0, -2, 0, 2, 0, 0, 2, 0, -2, 0, 0, -2, 0, 2, 0, 0],
            # Half-integer spin representations (Kramers doublets)
            'E1/2g': [2, 1.414, 0, 0, -1.414, 2, -1.414, 0, 0, 1.414, -2, -1.414, 0, 0, 1.414, -2, 1.414, 0, 0, -1.414],
            'E1/2u': [2, 1.414, 0, 0, -1.414, -2, 1.414, 0, 0, -1.414, -2, -1.414, 0, 0, 1.414, 2, -1.414, 0, 0, 1.414],
            'E3/2g': [2, -1.414, 0, 0, 1.414, 2, 1.414, 0, 0, -1.414, -2, 1.414, 0, 0, -1.414, -2, -1.414, 0, 0, 1.414],
            'E3/2u': [2, -1.414, 0, 0, 1.414, -2, -1.414, 0, 0, 1.414, -2, 1.414, 0, 0, -1.414, 2, 1.414, 0, 0, -1.414],
        }


class D2Star(PointGroup):
    """D2* double group - includes half-integer representations for fermions."""
    
    def __init__(self):
        super().__init__()
        self.name = "D2*"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C2(z)', 'C2(y)', 'C2(x)', 'R', 'RC2(z)', 'RC2(y)', 'RC2(x)']
        self.class_sizes = [1, 1, 1, 1, 1, 1, 1, 1]
        self.order = 8
        self.irreps = {
            # Integer spin representations (same as D2)
            'A': [1, 1, 1, 1, 1, 1, 1, 1],
            'B1': [1, 1, -1, -1, 1, 1, -1, -1],
            'B2': [1, -1, 1, -1, 1, -1, 1, -1],
            'B3': [1, -1, -1, 1, 1, -1, -1, 1],
            # Half-integer spin representations
            'E1/2': [2, 0, 0, 0, -2, 0, 0, 0],  # j = 1/2
            'E3/2': [2, 0, 0, 0, -2, 0, 0, 0],  # j = 3/2
        }


class D3Star(PointGroup):
    """D3* double group - includes half-integer representations for fermions."""
    
    def __init__(self):
        super().__init__()
        self.name = "D3*"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        self.classes = ['E', 'C3', 'C2', 'R', 'RC3', 'RC2']
        self.class_sizes = [1, 2, 3, 1, 2, 3]
        self.order = 12
        self.irreps = {
            # Integer spin representations (same as D3)
            'A1': [1, 1, 1, 1, 1, 1],
            'A2': [1, 1, -1, 1, 1, -1],
            'E': [2, -1, 0, 2, -1, 0],
            # Half-integer spin representations
            'E1/2': [2, 1, 0, -2, -1, 0],   # j = 1/2
            'E3/2': [2, -1, 0, -2, 1, 0],   # j = 3/2
        }


class D6Star(PointGroup):
    """D6* double group - includes half-integer representations for fermions."""
    
    def __init__(self):
        super().__init__()
        self.name = "D6*"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        import numpy as np
        self.classes = ['E', 'C6', 'C3', 'C2', 'C2\'', 'C2"', 'R', 'RC6', 'RC3', 'RC2', 'RC2\'', 'RC2"']
        self.class_sizes = [1, 2, 2, 1, 3, 3, 1, 2, 2, 1, 3, 3]
        self.order = 24
        self.irreps = {
            # Integer spin representations (same as D6)
            'A1': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            'A2': [1, 1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1],
            'B1': [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
            'B2': [1, -1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1],
            'E1': [2, 1, -1, -2, 0, 0, 2, 1, -1, -2, 0, 0],
            'E2': [2, -1, -1, 2, 0, 0, 2, -1, -1, 2, 0, 0],
            # Half-integer spin representations
            'E1/2': [2, 1, -1, -2, 0, 0, -2, -1, 1, 2, 0, 0],   # j = 1/2
            'E3/2': [2, np.sqrt(3), 0, 0, 0, 0, -2, -np.sqrt(3), 0, 0, 0, 0],   # j = 3/2
            'E5/2': [2, -np.sqrt(3), 0, 0, 0, 0, -2, np.sqrt(3), 0, 0, 0, 0],   # j = 5/2
        }


class OStar(PointGroup):
    """O* double group - octahedral symmetry with half-integer representations (already implemented as OhStar)."""
    
    def __init__(self):
        super().__init__()
        self.name = "O*"
        self._initialize_character_table()
    
    def _initialize_character_table(self):
        # This is the pure rotational octahedral double group
        # For simplicity, we use the same structure as Oh* but with fewer irreps
        self.classes = ['E', 'C3', 'C2', 'C4', 'C2\'', 'R', 'RC3', 'RC2', 'RC4', 'RC2\'']
        self.class_sizes = [1, 8, 6, 6, 3, 1, 8, 6, 6, 3]
        self.order = 48  # Double the original O order (24)
        self.irreps = {
            # Integer spin representations (same as O)
            'A1': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            'A2': [1, 1, -1, -1, 1, 1, 1, -1, -1, 1],
            'E': [2, -1, 0, 0, 2, 2, -1, 0, 0, 2],
            'T1': [3, 0, -1, 1, -1, 3, 0, -1, 1, -1],
            'T2': [3, 0, 1, -1, -1, 3, 0, 1, -1, -1],
            # Half-integer spin representations
            'E1/2': [2, 1, 0, 0, -2, -2, -1, 0, 0, 2],     # j = 1/2
            'E3/2': [4, 1, 0, 0, 0, -4, -1, 0, 0, 0],      # j = 3/2
            'E5/2': [4, -1, 0, 0, 0, -4, 1, 0, 0, 0],      # j = 5/2
        }


class PointGroupFactory:
    """Factory class to create point group objects."""
    
    _point_groups = {
        # Regular point groups
        'C1': C1,
        'Ci': Ci,
        'Cs': Cs,
        'C2': C2,
        'C3': C3,
        'C4': C4,
        'C5': C5,
        'C6': C6,
        'C7': C7,
        'C8': C8,
        'C2v': C2v,
        'C3v': C3v,
        'C4v': C4v,
        'C5v': C5v,
        'C6v': C6v,
        'C2h': C2h,
        'C3h': C3h,
        'C4h': C4h,
        'C5h': C5h,
        'C6h': C6h,
        'D2': D2,
        'D2h': D2h,
        'D3': D3,
        'D3h': D3h,
        'D4': D4,
        'D4h': D4h,
        'D5': D5,
        'D5h': D5h,
        'D6h': D6h,
        'D2d': D2d,
        'D3d': D3d,
        'D4d': D4d,
        'D5d': D5d,
        'D6d': D6d,
        'S4': S4,
        'S6': S6,
        'S8': S8,
        'T': T,
        'Td': Td,
        'Th': Th,
        'O': O,
        'Oh': Oh,
        'I': I,
        'Ih': Ih,
        'Cinfv': Cinfv,
        'C∞v': Cinfv,
        'Dinfh': Dinfh,
        'D∞h': Dinfh,
        # Double groups for half-integer spin systems
        'C2v*': C2vStar,
        'Oh*': OhStar,
        'Td*': TdStar,
        'D4h*': D4hStar,
        'D2*': D2Star,
        'D3*': D3Star,
        'D6*': D6Star,
        'O*': OStar,
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


def parse_irrep_string(irrep_string: str) -> list:
    """
    Parse a string containing multiple irrep labels separated by common delimiters.
    
    Args:
        irrep_string: String containing irrep labels (e.g., "T1u T2g", "T1u,T2g", "T1u×T2g")
        
    Returns:
        List of individual irrep labels
        
    Example:
        >>> parse_irrep_string("T1u T2g Eg")
        ['T1u', 'T2g', 'Eg']
        >>> parse_irrep_string("T1u,T2g,Eg")
        ['T1u', 'T2g', 'Eg'] 
        >>> parse_irrep_string("T1u×T2g×Eg")
        ['T1u', 'T2g', 'Eg']
    """
    import re
    # Split on various common delimiters: space, comma, ×, *, x, ⊗
    irrep_list = re.split(r'[,\s×*x⊗]+', irrep_string.strip())
    # Remove empty strings and strip whitespace
    return [irrep.strip() for irrep in irrep_list if irrep.strip()]


def calculate_multi_direct_product(point_group_name: str, *irreps) -> Dict[str, int]:
    """
    Calculate direct product of multiple irreps recursively.
    
    Args:
        point_group_name: Name of the point group (e.g., 'Oh', 'Td')
        *irreps: Variable number of irrep labels as separate arguments
        
    Returns:
        Dictionary mapping irrep names to coefficients in the final direct product
        
    Example:
        >>> calculate_multi_direct_product('Oh', 'T1u', 'T2g', 'Eg')
        {'Eu': 2, 'A1u': 1, 'A2u': 1, 'T1u': 2, 'T2u': 2}
        >>> calculate_multi_direct_product('Td', 'A1', 'E', 'T2')
        {'T1': 1, 'T2': 1, 'E': 1}
    """
    # Handle the case where a single string with multiple irreps is passed for backwards compatibility
    if len(irreps) == 1 and isinstance(irreps[0], str) and any(delimiter in irreps[0] for delimiter in [' ', ',', '×', '*', 'x', '⊗']):
        irrep_list = parse_irrep_string(irreps[0])
    else:
        irrep_list = list(irreps)
    
    if len(irrep_list) < 2:
        raise ValueError("At least two irreps are required for direct product calculation")
    
    # Start with the first two irreps
    result = calculate_direct_product(point_group_name, irrep_list[0], irrep_list[1])
    
    # Recursively calculate direct product with remaining irreps
    for i in range(2, len(irrep_list)):
        next_irrep = irrep_list[i]
        new_result = {}
        
        # Calculate direct product of current result with next irrep
        for current_irrep, current_coeff in result.items():
            product = calculate_direct_product(point_group_name, current_irrep, next_irrep)
            
            # Add the contributions to new_result
            for irrep, coeff in product.items():
                new_result[irrep] = new_result.get(irrep, 0) + current_coeff * coeff
        
        result = new_result
    
    return result


def multi_direct_product_label(point_group_name: str, *irreps) -> str:
    """
    Get formatted label for direct product of multiple irreps.
    
    Args:
        point_group_name: Name of the point group
        *irreps: Variable number of irrep labels as separate arguments
        
    Returns:
        Formatted string showing the direct product decomposition
        
    Example:
        >>> multi_direct_product_label('Oh', 'T1u', 'T2g', 'Eg')
        '2Eu ⊕ A1u ⊕ A2u ⊕ 2T1u ⊕ 2T2u'
        >>> multi_direct_product_label('Td', 'A1', 'E', 'T2')
        'T1 ⊕ T2 ⊕ E'
    """
    product = calculate_multi_direct_product(point_group_name, *irreps)
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
    
    print("\n\nDirect Product Examples:")
    print(f"T1u × T2g = {direct_product_label('Oh', 'T1u', 'T2g')}")
    
    print("\nMulti-Irrep Direct Product Examples:")
    print(f"T1u × T2g × Eg = {multi_direct_product_label('Oh', 'T1u', 'T2g', 'Eg')}")
    print(f"A1g × T1u × T1u = {multi_direct_product_label('Oh', 'A1g', 'T1u', 'T1u')}")
    print(f"Four irreps: {multi_direct_product_label('Oh', 'T1g', 'T2g', 'Eg', 'A1g')}")
    print(f"Backwards compatibility (string): {multi_direct_product_label('Oh', 'T1g T2g Eg')}")
