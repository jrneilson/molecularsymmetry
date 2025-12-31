#!/usr/bin/env python3
"""
Focused Demo: Octahedral Complex Symmetry Analysis

This demonstrates the molecular_symmetry package with a detailed
octahedral ML6 complex example.
"""

from molecular_symmetry import get_point_group

def main():
    print("=" * 80)
    print("Molecular Orbital Symmetry Analysis: Octahedral ML6 Complex")
    print("=" * 80)
    
    # Get the Oh point group
    oh = get_point_group('Oh')
    
    # Display the character table
    oh.print_character_table()
    
    print("\n\n" + "=" * 80)
    print("REDUCIBLE REPRESENTATIONS FOR OCTAHEDRAL COMPLEX")
    print("=" * 80)
    
    # Define various orbital sets and their representations
    orbital_sets = {
        "6 Ligand σ orbitals": {
            "rep": [6, 0, 0, 2, 2, 0, 0, 0, 4, 2],
            "explanation": "Six sigma-donor orbitals from the ligands"
        },
        "12 Ligand π orbitals": {
            "rep": [12, 0, 0, 0, -4, 0, 0, 0, 0, 0],
            "explanation": "Twelve pi orbitals (2 per ligand, perpendicular to M-L bond)"
        },
        "Metal s orbital": {
            "rep": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            "explanation": "Single metal s orbital"
        },
        "Metal p orbitals": {
            "rep": [3, 0, -1, 1, -1, -3, -1, 0, 1, 1],
            "explanation": "Three metal p orbitals (px, py, pz)"
        },
        "Metal d orbitals": {
            "rep": [5, -1, 1, -1, 1, 5, -1, -1, 1, 1],
            "explanation": "Five metal d orbitals"
        },
        "Metal f orbitals": {
            "rep": [7, 1, -1, -1, -1, -7, 1, -1, 1, 1],
            "explanation": "Seven metal f orbitals (rarely involved in bonding)"
        }
    }
    
    # Analyze each orbital set
    for name, data in orbital_sets.items():
        print(f"\n{name}:")
        print(f"  Description: {data['explanation']}")
        print(f"  Γ = {data['rep']}")
        
        # Reduce the representation
        irreps = oh.reduce_representation(data['rep'])
        label = oh.get_symmetry_label(data['rep'])
        
        print(f"  Irreducible representations: {irreps}")
        print(f"  Symmetry label: {label}")
    
    # Detailed analysis
    print("\n\n" + "=" * 80)
    print("MOLECULAR ORBITAL DIAGRAM INTERPRETATION")
    print("=" * 80)
    
    print("\n1. LIGAND ORBITALS:")
    print("-" * 80)
    print("   Sigma bonding (A1g + Eg + T1u):")
    print("   • A1g: totally symmetric combination")
    print("   • Eg: doubly degenerate combination")
    print("   • T1u: triply degenerate combination")
    print()
    print("   Pi bonding (T1g + T2g + T1u + T2u):")
    print("   • T1g, T2g: combinations with gerade (g) symmetry")
    print("   • T1u, T2u: combinations with ungerade (u) symmetry")
    
    print("\n2. METAL ORBITALS:")
    print("-" * 80)
    print("   s orbital → A1g (spherically symmetric)")
    print("   p orbitals → T1u (can mix with ligand T1u for bonding)")
    print("   d orbitals → Eg + T2g (crystal field splitting!)")
    print("     • Eg: dz² and dx²-y² (point toward ligands, σ-antibonding)")
    print("     • T2g: dxy, dxz, dyz (point between ligands, π-bonding)")
    
    print("\n3. MOLECULAR ORBITAL FORMATION:")
    print("-" * 80)
    print("   σ-bonding MOs:")
    print("   • Ligand A1g + Metal A1g (s) → bonding & antibonding A1g")
    print("   • Ligand Eg + Metal Eg (dz², dx²-y²) → bonding & antibonding Eg")
    print("   • Ligand T1u + Metal T1u (p orbitals) → bonding & antibonding T1u")
    print()
    print("   π-bonding MOs (if ligands are π-donors):")
    print("   • Ligand T2g + Metal T2g (dxy, dxz, dyz) → affects Δo")
    print()
    print("   Result: The famous crystal field splitting!")
    print("   • eg* orbitals (higher energy): dz², dx²-y²")
    print("   • t2g orbitals (lower energy): dxy, dxz, dyz")
    print("   • Δo (octahedral splitting) = E(eg*) - E(t2g)")
    
    print("\n4. SYMMETRY SELECTION RULES:")
    print("-" * 80)
    print("   Orbitals can only interact if they have the same symmetry!")
    print("   • Metal A1g (s) ↔ Ligand A1g σ-combinations")
    print("   • Metal Eg (dz², dx²-y²) ↔ Ligand Eg σ-combinations")
    print("   • Metal T1u (px, py, pz) ↔ Ligand T1u σ-combinations")
    print("   • Metal T2g (dxy, dxz, dyz) ↔ Ligand T2g π-combinations")
    
    print("\n\n" + "=" * 80)
    print("QUICK REFERENCE")
    print("=" * 80)
    print("\nHow to use this package:")
    print("  from molecular_symmetry import get_point_group")
    print("  oh = get_point_group('Oh')")
    print("  irreps = oh.reduce_representation([6, 0, 0, 2, 2, 0, 0, 0, 4, 2])")
    print("  print(oh.get_symmetry_label([6, 0, 0, 2, 2, 0, 0, 0, 4, 2]))")
    print("\nAvailable point groups:")
    from molecular_symmetry import PointGroupFactory
    groups = PointGroupFactory.list_available()
    print(f"  {', '.join(groups[:15])}")
    print(f"  ... and {len(groups) - 15} more!")
    
    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
