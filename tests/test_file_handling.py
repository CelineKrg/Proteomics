from proteomics.file_handling import read_fasta
def test_read_fasta():
    tmp_fasta_path = 'C:\\Users\\Celine\Desktop\\University\\S5\\advanced-python_course\\data\\dummy-proteins.fasta'
    protein_map = read_fasta(tmp_fasta_path)
    
    # Replace the strings with your fasta content
    # which you expect to be now available as a dictionary
    assert protein_map["34G7D2"] == "DWRGDSJKLSDFJ"
    assert protein_map["D49J2K"] == "KFPASOKFVVDSKFSL"

test_read_fasta()