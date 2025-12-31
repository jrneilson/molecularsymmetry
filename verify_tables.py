#!/usr/bin/env python3
"""
Verify all character tables for correctness using orthogonality relations.
"""

import numpy as np
from molecular_symmetry import PointGroupFactory

def verify_character_table(pg):
    """Verify a character table using orthogonality relations."""
    print(f"\n{'='*70}")
    print(f"Verifying {pg.name} character table")
    print(f"{'='*70}")
    
    # Check 1: Sum of squares for each irrep should equal the group order
    print("\nOrthonormality check (Σ n·χ² = h):")
    all_good = True
    for name, chars in pg.irreps.items():
        sum_sq = sum(n * abs(c)**2 for n, c in zip(pg.class_sizes, chars))
        status = "✓" if abs(sum_sq - pg.order) < 0.001 else "✗"
        if status == "✗":
            all_good = False
        print(f"  {name:8s}: {sum_sq:6.1f} (expected {pg.order}) {status}")
    
    # Check 2: Orthogonality between different irreps
    print("\nOrthogonality check (different irreps):")
    names = list(pg.irreps.keys())
    for i, name1 in enumerate(names):
        for name2 in names[i+1:]:
            dot = sum(n * c1 * np.conj(c2) for n, c1, c2 in 
                     zip(pg.class_sizes, pg.irreps[name1], pg.irreps[name2]))
            if abs(dot) > 0.001:
                print(f"  {name1} · {name2} = {dot.real:.3f} (should be 0!) ✗")
                all_good = False
    
    if all_good:
        print("  All pairs orthogonal ✓")
    
    # Check 3: Number of irreps
    n_irreps = len(pg.irreps)
    n_classes = len(pg.classes)
    status = "✓" if n_irreps == n_classes else "✗"
    print(f"\nNumber of irreps: {n_irreps}, Number of classes: {n_classes} {status}")
    
    # Check 4: Sum of squared dimensions
    dim_sum = sum((chars[0])**2 for chars in pg.irreps.values())
    status = "✓" if abs(dim_sum - pg.order) < 0.001 else "✗"
    print(f"Σ(dimension²): {dim_sum:.1f} (expected {pg.order}) {status}")
    
    return all_good


def main():
    print("CHARACTER TABLE VERIFICATION SUITE")
    print("="*70)
    
    # Get all point groups
    pg_names = PointGroupFactory.list_available()
    
    # Skip infinite groups
    skip_groups = ['Cinfv', 'C∞v', 'Dinfh', 'D∞h']
    pg_names = [name for name in pg_names if name not in skip_groups]
    
    results = {}
    for name in pg_names:
        try:
            pg = PointGroupFactory.create(name)
            is_valid = verify_character_table(pg)
            results[name] = is_valid
        except Exception as e:
            print(f"\n{'='*70}")
            print(f"ERROR in {name}: {e}")
            print(f"{'='*70}")
            results[name] = False
    
    # Summary
    print("\n\n" + "="*70)
    print("VERIFICATION SUMMARY")
    print("="*70)
    
    passed = [name for name, valid in results.items() if valid]
    failed = [name for name, valid in results.items() if not valid]
    
    print(f"\nPassed: {len(passed)}/{len(results)}")
    if failed:
        print(f"Failed: {', '.join(failed)}")
    else:
        print("All character tables are mathematically correct! ✓")


if __name__ == "__main__":
    main()
