import re
from .cleavage_patterns import enzyme_cleavage_patterns

def digest_protein_sequence(protein_seq, cleave_pattern):
    """
    function to simulate digestion of protein sequences 

    Parameters
    ----------
    protein_seq : str
        full amino-acid sequence of the protein to be digested
    cleave_pattern : str
        regular-expression pattern defining the positions where the sequence will be split
   
    Returns
    -------
    peptides : list of str
        list of peptide fragments produced after applying the cleavage pattern to the protein sequence

    """
    peptides = re.split(cleave_pattern, protein_seq)

    return peptides, f'Nr. of digested peptides: {len(peptides)}'


def digest_protein_collection(protein_map, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    function to simulate digestion of proteins. 

    Parameters
    ----------
    protein_map : dict
        dictionary mapping protein ID to the corresponding amino-acid sequences
    cleave_pattern : str or Pattern
        regular-expression pattern defining the positions where the sequence will be split
    min_pep_len : int, optional
        Minimum allowed peptide length after digestion with the default set to 5
    max_pep_len : int, optional
        Maximum allowed peptide length after digestion with the default set to 30
        
       
    Returns
    -------
   digested_proteins : dict
        dictionary containing the peptide fragments produced after applying the cleavage pattern to the protein sequence with the length in the defined range and the number of retained peptides
    """
    if min_pep_len > max_pep_len: 
        min_pep_len, max_pep_len = max_pep_len, min_pep_len

    digested_proteins = {}


    for protein_id, protein_seq in protein_map.items():
        if not protein_seq: 
            digested_proteins[protein_id] = ([], 'Nr. of digested peptides: 0')
            continue
        
        peptides = re.split(cleave_pattern, protein_seq)
        filtered = [p for p in peptides if min_pep_len <= len(p) <= max_pep_len] 
        info = f'Nr. of digested peptides: {len(filtered)}'
        digested_proteins[protein_id] = (filtered, info)

    return digested_proteins


def compute_sequence_coverage(protein_seq, peptides):
    """
    function to compute the sequence coverage of a protein with its detected peptides

    Parameters
    protein_seq : str
        list containing the complete protein_sequence before the digestion 
    peptides : list
        list of peptide fragments produced after applying the cleavage pattern to the protein sequence with the length in the defined range
    ----------
    Returns
    coverage_percent : float
        protein coverage (in percent)
    -------
    """
    if not protein_seq or not peptides: # if one of them is empty = 0 
        return 0
    coverage = set() # set to save each element only one time -> no doubles
    for pep in peptides: 
        idx_start = 0 # start searching for the peptide at the start of proteinsequence 
        while True: # because we dont know how often the peptide is in the sequence
            idx = protein_seq.find(pep, idx_start) # gives the first position where the peptide is found
            if idx == -1: # if the peptide is not found: leave the loop
                break 
            coverage.update(range(idx, idx + len(pep)))  # add all indexes of the found peptide to the set (coverage contains all indexes of the protein_seq that are covered by the peptides)
            idx_start = idx_start + 1 # to cover all indices of the peptide (if it is multiple times in the protein sequence)
    coverage_percent = len(coverage) / len(protein_seq) * 100 # get the percentage (all the indices that are in the protein_seq that are covered by peptides / complete length of preotein_seq)

    return coverage_percent