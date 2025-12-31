#!/usr/bin/env python3
"""
Example usage of the molecular symmetry package.

This script demonstrates how to use the package for various molecular systems.
"""

from molecular_symmetry import get_point_group, PointGroupFactory


def example_octahedral_complex():
    """Example: ML6 octahedral complex."""
    print("=" * 80)
    print("EXAMPLE 1: Octahedral ML6 Complex (Oh symmetry)")
    print("=" * 80)
    
    oh = get_point_group('Oh')
    oh.print_character_table()
    
    examples = {
        "Sigma bonding (6 ligand σ orbitals)": [6, 0, 0, 2, 2, 0, 0, 0, 4, 2],
        "Pi bonding (12 ligand π orbitals)": [12, 0, 0, 0, -4, 0, 0, 0, 0, 0],
        "Metal s orbital": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        "Metal p orbitals": [3, 0, -1, 1, -1, -3, -1, 0, 1, 1],
        "Metal d orbitals": [5, -1, 1, -1, 1, 5, -1, -1, 1, 1],
        "Metal f orbitals": [7, 1, -1, -1, -1, -7, 1, -1, 1, 1],
    }
    
    print("\n\nMolecular Orbital Analysis:")
    print("-" * 80)
    for description, rep in examples.items():
        irreps = oh.reduce_representation(rep)
        label = oh.get_symmetry_label(rep)
        print(f"\n{description}:")
        print(f"  Γ = {rep}")
        print(f"  → {label}")


def example_tetrahedral_complex():
    """Example: ML4 tetrahedral complex."""
    print("\n\n" + "=" * 80)
    print("EXAMPLE 2: Tetrahedral ML4 Complex (Td symmetry)")
    print("=" * 80)
    
    td = get_point_group('Td')
    td.print_character_table()
    
    examples = {
        "Sigma bonding (4 ligand σ orbitals)": [4, 1, 0, 0, 2],
        "Metal s orbital": [1, 1, 1, 1, 1],
        "Metal p orbitals": [3, 0, -1, -1, 1],
        "Metal d orbitals": [5, -1, 1, -1, 1],
    }
    
    print("\n\nMolecular Orbital Analysis:")
    print("-" * 80)
    for description, rep in examples.items():
        label = td.get_symmetry_label(rep)
        print(f"\n{description}:")
        print(f"  → {label}")


def example_square_planar_complex():
    """Example: ML4 square planar complex."""
    print("\n\n" + "=" * 80)
    print("EXAMPLE 3: Square Planar ML4 Complex (D4h symmetry)")
    print("=" * 80)
    
    d4h = get_point_group('D4h')
    d4h.print_character_table()
    
    examples = {
        "Sigma bonding (4 ligand σ orbitals)": [4, 0, 0, 2, 0, 0, 0, 4, 2, 0],
        "Metal s orbital": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        "Metal px, py orbitals": [2, 0, 0, 0, -2, -2, 0, 0, 0, 2],
        "Metal pz orbital": [1, 1, 1, -1, -1, -1, -1, -1, 1, 1],
        "Metal dxy, dx²-y² orbitals": [2, 0, 2, 0, -2, 2, 0, 2, 2, -2],
        "Metal dxz, dyz orbitals": [2, 0, 0, 0, 2, -2, 0, 0, 0, -2],
        "Metal dz² orbital": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    }
    
    print("\n\nMolecular Orbital Analysis:")
    print("-" * 80)
    for description, rep in examples.items():
        label = d4h.get_symmetry_label(rep)
        print(f"\n{description}:")
        print(f"  → {label}")


def example_water_molecule():
    """Example: H2O water molecule."""
    print("\n\n" + "=" * 80)
    print("EXAMPLE 4: Water Molecule H2O (C2v symmetry)")
    print("=" * 80)
    
    c2v = get_point_group('C2v')
    c2v.print_character_table()
    
    examples = {
        "Two H 1s orbitals": [2, 0, 2, 0],
        "O 2s orbital": [1, 1, 1, 1],
        "O 2px orbital": [1, -1, 1, -1],
        "O 2py orbital": [1, -1, -1, 1],
        "O 2pz orbital": [1, 1, 1, 1],
    }
    
    print("\n\nMolecular Orbital Analysis:")
    print("-" * 80)
    for description, rep in examples.items():
        label = c2v.get_symmetry_label(rep)
        print(f"\n{description}:")
        print(f"  → {label}")


def example_ammonia_molecule():
    """Example: NH3 ammonia molecule."""
    print("\n\n" + "=" * 80)
    print("EXAMPLE 5: Ammonia Molecule NH3 (C3v symmetry)")
    print("=" * 80)
    
    c3v = get_point_group('C3v')
    c3v.print_character_table()
    
    examples = {
        "Three H 1s orbitals": [3, 0, 1],
        "N 2s orbital": [1, 1, 1],
        "N 2px, 2py orbitals": [2, -1, 0],
        "N 2pz orbital": [1, 1, 1],
    }
    
    print("\n\nMolecular Orbital Analysis:")
    print("-" * 80)
    for description, rep in examples.items():
        label = c3v.get_symmetry_label(rep)
        print(f"\n{description}:")
        print(f"  → {label}")


def example_benzene_molecule():
    """Example: C6H6 benzene molecule."""
    print("\n\n" + "=" * 80)
    print("EXAMPLE 6: Benzene Molecule C6H6 (D6h symmetry)")
    print("=" * 80)
    
    d6h = get_point_group('D6h')
    d6h.print_character_table()
    
    examples = {
        "Six C pz orbitals (π system)": [6, 0, 0, 0, 2, 0, 0, 0, 0, 6, 0, 2],
        "Six C-H σ bonds": [6, 0, 0, 0, 2, 4, 0, 0, 0, 6, 4, 2],
    }
    
    print("\n\nMolecular Orbital Analysis:")
    print("-" * 80)
    for description, rep in examples.items():
        label = d6h.get_symmetry_label(rep)
        print(f"\n{description}:")
        print(f"  → {label}")


def example_methane_molecule():
    """Example: CH4 methane molecule."""
    print("\n\n" + "=" * 80)
    print("EXAMPLE 7: Methane Molecule CH4 (Td symmetry)")
    print("=" * 80)
    
    td = get_point_group('Td')
    td.print_character_table()
    
    examples = {
        "Four H 1s orbitals": [4, 1, 0, 0, 2],
        "C 2s orbital": [1, 1, 1, 1, 1],
        "C 2p orbitals": [3, 0, -1, -1, 1],
    }
    
    print("\n\nMolecular Orbital Analysis:")
    print("-" * 80)
    for description, rep in examples.items():
        label = td.get_symmetry_label(rep)
        print(f"\n{description}:")
        print(f"  → {label}")
    
    print("\n\nInterpretation:")
    print("-" * 80)
    print("• Four H 1s orbitals → A1 ⊕ T2")
    print("• C 2s orbital → A1")
    print("• C 2p orbitals → T2")
    print("• SALCs for bonding: A1 (from H orbitals) + A1 (C 2s)")
    print("                      T2 (from H orbitals) + T2 (C 2p)")
    print("• This gives four equivalent sp³ hybrid orbitals!")


def example_ferrocene():
    """Example: Fe(C5H5)2 ferrocene (staggered)."""
    print("\n\n" + "=" * 80)
    print("EXAMPLE 8: Ferrocene Fe(C5H5)2 - Staggered (D5d symmetry)")
    print("=" * 80)
    
    d5h = get_point_group('D5h')
    d5h.print_character_table()
    
    print("\n\nNote: Ferrocene in eclipsed conformation has D5h symmetry")
    print("      Staggered conformation has D5d symmetry")
    print("\nFor D5h (eclipsed):")
    
    examples = {
        "Fe d orbitals (simplified)": [5, 0, 0, -1, 5, 0, 0, -1],
        "Two Cp π systems (10 π MOs total)": [10, 0, 0, 0, 10, 0, 0, 0],
    }
    
    print("\nMolecular Orbital Analysis:")
    print("-" * 80)
    for description, rep in examples.items():
        label = d5h.get_symmetry_label(rep)
        print(f"\n{description}:")
        print(f"  → {label}")


def example_linear_molecules():
    """Examples of linear molecules."""
    print("\n\n" + "=" * 80)
    print("EXAMPLE 9: Linear Molecules")
    print("=" * 80)
    
    # CO molecule (C∞v)
    print("\nCO Molecule (C∞v symmetry - heteronuclear):")
    print("-" * 80)
    cinfv = get_point_group('C∞v')
    cinfv.print_character_table()
    
    # CO2 molecule (D∞h)
    print("\n\nCO2 Molecule (D∞h symmetry - homonuclear):")
    print("-" * 80)
    dinfh = get_point_group('D∞h')
    dinfh.print_character_table()
    
    print("\n\nNote: Linear molecules have infinite symmetry groups.")
    print("σ orbitals → Σ+ symmetry")
    print("π orbitals → Π symmetry")
    print("δ orbitals → Δ symmetry")


def list_all_point_groups():
    """List all available point groups."""
    print("\n\n" + "=" * 80)
    print("ALL AVAILABLE POINT GROUPS")
    print("=" * 80)
    
    groups = PointGroupFactory.list_available()
    
    categories = {
        "Simple groups": ['C1', 'Ci', 'Cs'],
        "Cyclic groups (Cn)": ['C2', 'C3'],
        "Pyramidal groups (Cnv)": ['C2v', 'C3v', 'C4v', 'C5v', 'C6v'],
        "Dihedral groups (Dn)": ['D2', 'D3'],
        "Dihedral with horizontal plane (Dnh)": ['D2h', 'D3h', 'D4h', 'D5h', 'D6h'],
        "Cubic groups": ['Td', 'Oh', 'Ih'],
        "Linear groups": ['C∞v', 'Cinfv', 'D∞h', 'Dinfh']
    }
    
    for category, group_list in categories.items():
        available = [g for g in group_list if g in groups]
        if available:
            print(f"\n{category}:")
            print(f"  {', '.join(available)}")


if __name__ == "__main__":
    # Run all examples
    example_octahedral_complex()
    example_tetrahedral_complex()
    example_square_planar_complex()
    example_water_molecule()
    example_ammonia_molecule()
    example_benzene_molecule()
    example_methane_molecule()
    example_ferrocene()
    example_linear_molecules()
    list_all_point_groups()
    
    print("\n\n" + "=" * 80)
    print("Examples complete!")
    print("=" * 80)
    print("\nTo use in your own code:")
    print("  from molecular_symmetry import get_point_group")
    print("  pg = get_point_group('Oh')")
    print("  irreps = pg.reduce_representation([6, 0, 0, 2, 2, 0, 0, 0, 4, 2])")
    print("  print(pg.get_symmetry_label([6, 0, 0, 2, 2, 0, 0, 0, 4, 2]))")
