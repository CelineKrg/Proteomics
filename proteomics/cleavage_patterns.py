enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)', # digest after Lysin
    'LysN': r'(?=K)', # digest before Lysin
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)', # [] means this or that
}