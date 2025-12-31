#!/usr/bin/env python3
"""
Test all orbital representations in examples for correctness.
"""

from molecular_symmetry import get_point_group

def test_oh_orbitals():
    """Test Oh (octahedral) orbital representations."""
    print("\n" + "="*70)
    print("Testing Oh (Octahedral) Orbitals")
    print("="*70)
    
    oh = get_point_group('Oh')
    
    tests = {
        "6 σ orbitals": {
            "rep": [6, 0, 0, 2, 2, 0, 0, 0, 4, 2],
            "expected": "A1g ⊕ Eg ⊕ T1u"
        },
        "12 π orbitals": {
            "rep": [12, 0, 0, 0, -4, 0, 0, 0, 0, 0],
            "expected": "T1g ⊕ T2g ⊕ T1u ⊕ T2u"
        },
        "Metal s": {
            "rep": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            "expected": "A1g"
        },
        "Metal p": {
            "rep": [3, 0, -1, 1, -1, -3, -1, 0, 1, 1],
            "expected": "T1u"
        },
        "Metal d": {
            "rep": [5, -1, 1, -1, 1, 5, -1, -1, 1, 1],
            "expected": "Eg ⊕ T2g"
        },
        "Metal f": {
            "rep": [7, 1, -1, -1, -1, -7, 1, -1, 1, 1],
            "expected": "A2u ⊕ T1u ⊕ T2u"
        },
    }
    
    all_pass = True
    for name, data in tests.items():
        result = oh.get_symmetry_label(data['rep'])
        status = "✓" if result == data['expected'] else "✗"
        if status == "✗":
            all_pass = False
        print(f"{name:15s}: {result:30s} (expected: {data['expected']}) {status}")
    
    return all_pass


def test_td_orbitals():
    """Test Td (tetrahedral) orbital representations."""
    print("\n" + "="*70)
    print("Testing Td (Tetrahedral) Orbitals")
    print("="*70)
    
    td = get_point_group('Td')
    
    tests = {
        "4 σ orbitals": {
            "rep": [4, 1, 0, 0, 2],
            "expected": "A1 ⊕ T2"
        },
        "Metal s": {
            "rep": [1, 1, 1, 1, 1],
            "expected": "A1"
        },
        "Metal p": {
            "rep": [3, 0, -1, -1, 1],
            "expected": "T2"
        },
        "Metal d": {
            "rep": [5, -1, 1, -1, 1],
            "expected": "E ⊕ T2"
        },
    }
    
    all_pass = True
    for name, data in tests.items():
        result = td.get_symmetry_label(data['rep'])
        status = "✓" if result == data['expected'] else "✗"
        if status == "✗":
            all_pass = False
        print(f"{name:15s}: {result:30s} (expected: {data['expected']}) {status}")
    
    return all_pass


def test_d4h_orbitals():
    """Test D4h (square planar) orbital representations."""
    print("\n" + "="*70)
    print("Testing D4h (Square Planar) Orbitals")
    print("="*70)
    
    d4h = get_point_group('D4h')
    
    tests = {
        "4 σ orbitals": {
            "rep": [4, 0, 0, 2, 0, 0, 0, 4, 2, 0],
            "expected": "A1g ⊕ B1g ⊕ Eu"
        },
        "Metal s": {
            "rep": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            "expected": "A1g"
        },
        "Metal d": {
            "rep": [5, -1, 1, 1, 1, 5, -1, 1, 1, 1],
            "expected": "A1g ⊕ B1g ⊕ B2g ⊕ Eg"
        },
    }
    
    all_pass = True
    for name, data in tests.items():
        result = d4h.get_symmetry_label(data['rep'])
        status = "✓" if result == data['expected'] else "✗"
        if status == "✗":
            all_pass = False
        print(f"{name:15s}: {result:30s} (expected: {data['expected']}) {status}")
    
    return all_pass


def test_c2v_orbitals():
    """Test C2v (water) orbital representations."""
    print("\n" + "="*70)
    print("Testing C2v (Water) Orbitals")
    print("="*70)
    
    c2v = get_point_group('C2v')
    
    tests = {
        "2 H 1s": {
            "rep": [2, 0, 2, 0],
            "expected": "A1 ⊕ B1"
        },
        "O 2s": {
            "rep": [1, 1, 1, 1],
            "expected": "A1"
        },
        "O 2px": {
            "rep": [1, -1, 1, -1],
            "expected": "B1"
        },
        "O 2py": {
            "rep": [1, -1, -1, 1],
            "expected": "B2"
        },
        "O 2pz": {
            "rep": [1, 1, 1, 1],
            "expected": "A1"
        },
    }
    
    all_pass = True
    for name, data in tests.items():
        result = c2v.get_symmetry_label(data['rep'])
        status = "✓" if result == data['expected'] else "✗"
        if status == "✗":
            all_pass = False
        print(f"{name:15s}: {result:30s} (expected: {data['expected']}) {status}")
    
    return all_pass


def test_c3v_orbitals():
    """Test C3v (ammonia) orbital representations."""
    print("\n" + "="*70)
    print("Testing C3v (Ammonia) Orbitals")
    print("="*70)
    
    c3v = get_point_group('C3v')
    
    tests = {
        "3 H 1s": {
            "rep": [3, 0, 1],
            "expected": "A1 ⊕ E"
        },
        "N 2s": {
            "rep": [1, 1, 1],
            "expected": "A1"
        },
        "N 2px,2py": {
            "rep": [2, -1, 0],
            "expected": "E"
        },
        "N 2pz": {
            "rep": [1, 1, 1],
            "expected": "A1"
        },
    }
    
    all_pass = True
    for name, data in tests.items():
        result = c3v.get_symmetry_label(data['rep'])
        status = "✓" if result == data['expected'] else "✗"
        if status == "✗":
            all_pass = False
        print(f"{name:15s}: {result:30s} (expected: {data['expected']}) {status}")
    
    return all_pass


def main():
    print("="*70)
    print("ORBITAL REPRESENTATION VERIFICATION")
    print("="*70)
    
    results = []
    results.append(("Oh", test_oh_orbitals()))
    results.append(("Td", test_td_orbitals()))
    results.append(("D4h", test_d4h_orbitals()))
    results.append(("C2v", test_c2v_orbitals()))
    results.append(("C3v", test_c3v_orbitals()))
    
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    for name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{name:10s}: {status}")
    
    all_pass = all(passed for _, passed in results)
    if all_pass:
        print("\n✓ All orbital representations verified!")
    else:
        print("\n✗ Some orbital representations need fixing")
    
    return all_pass


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
