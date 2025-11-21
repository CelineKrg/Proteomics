from proteomics.mass_spectra_simulation import (
    calculate_mol_mass,
    calculate_mol_mass_collection, 
    calculate_mz_collection, 
    fragment_peptide
)                                
from proteomics.amino_acid_mass_dalton import amino_acid_mass_dalton

def test_calculate_mol_mass():
    aa_mass_dict = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}
    
    test_peptide1 = "ARGPEPKL"
    test_peptide2 = "CDFZERA"
    expected_mass1 = 849.01
    calculated_mass1 = calculate_mol_mass(test_peptide1, aa_mass_dict)
    
    assert calculated_mass1[test_peptide1] == expected_mass1

    try:
        calculate_mol_mass(test_peptide2, aa_mass_dict)
        assert False
    except ValueError:
        pass

test_calculate_mol_mass()


def test_calculate_mol_mass_collection():
    aa_mass_dict = amino_acid_mass_dalton
    peptides1 = ['ADFRK', 'PRKKL', 'MLDTI']
    peptides2 = ['ADFRK', 'PRKKL', 'MLDTI', 'WAZTC']
    expected1 = {
        'ADFRK': 617.71,
        'PRKKL': 622.81,
        'MLDTI': 573.71
    }

    actual1 = calculate_mol_mass_collection(peptides1, amino_acid_mass_dalton)

    assert actual1 == expected1
    
    try:
        calculate_mol_mass_collection(peptides2, amino_acid_mass_dalton)
        assert False
    except ValueError:
        pass

test_calculate_mol_mass_collection()


def test_calculate_mz_collection():
    peptide_mass_map1 = {
        'ADFRK': 617.71,
        'PRKKL': 622.81,
        'MLDTI': 573.71
    }

    peptide_mass_map2 = {
        'ADFRK': 617.71,
        'PRKKL': -23.1,
        'MLDTI': 573.71
    }

    actual1 = calculate_mz_collection(peptide_mass_map1, charge=2)
    expected1 = {
        'ADFRK': 309.862, 
        'PRKKL': 312.412, 
        'MLDTI': 287.862
        }

    assert actual1 == expected1

    try:
        calculate_mz_collection(peptide_mass_map1, charge=0)
        assert False
    except ValueError:
        pass

    try:
        calculate_mz_collection(peptide_mass_map2, charge=2)
        assert False
    except ValueError:
        pass

test_calculate_mz_collection()


def test_fragment_peptide():
    peptide = "ABCD"
    expected = ["A", "AB", "ABC", "BCD", "CD", "D"]
    actual = fragment_peptide(peptide)

    assert set(actual) == set(expected)

test_fragment_peptide()