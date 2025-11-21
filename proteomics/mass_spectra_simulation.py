import numpy as np
import matplotlib.pyplot as plt


def calculate_mol_mass(peptide_seq, amino_acid_mass_dict):
    """
    function to calculate the molecular mass of a single peptide sequence 

    Parameters
    peptide_seq = str 
        the peptide sequences
    
    amino_acid_mass_dict : dict
        dictionary containing the amino-acid one-letter codes (keys) and their masses in Daltons (values)
    ----------

    Returns
    peptide_mass : dict
        dictionary containing the peptides (keys) and their mass (values)
    -------

    """
    peptide_mass = {}
    mass = 0.0
    for aminoacid in peptide_seq: 
        if aminoacid not in amino_acid_mass_dict: 
            raise ValueError (f"Aminoacid not known: {aminoacid}")
        else: 
            mass += amino_acid_mass_dict[aminoacid]
    peptide_mass[peptide_seq] = mass
    return peptide_mass


def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
    Add a short description here.
    
    Parameters
    peptides = list  
        list containing the peptide sequences
    
    amino_acid_mass_dict : dict
        dictionary containing the amino-acid one-letter codes (keys) and their masses in Daltons (values)
    ----------

    Returns
    peptide_mass : dict
        dictionary containing the peptides (keys) and their mass (values)
    ----------

    Returns
    -------

    """

    peptide_mass = {}

    for pep in peptides: 
        mass = 0.0
        for aminoacid in pep: 
            if aminoacid not in amino_acid_mass_dict: 
                raise ValueError (f"Aminoacid not known: {aminoacid}")
            else: 
                mass += amino_acid_mass_dict[aminoacid]
        peptide_mass[pep] = mass
    return peptide_mass


def calculate_mz_collection(peptide_mass_map, charge=2, proton_mass=1.007):
    """
    function to map the peptides to their m/z value

    Parameters
    ----------
    peptide_mass_map : dict
        dictionary containing the peptides (keys) and their mass in Daltons (values)
    
    charge : int, optimal
        charge for the MS1 peptide with the default set to +2

    proton_mass : float, optimal 
        proton mass with the default set to 1.007

    Returns
    -------
    peptides_mz_value : dict
        dictionary containing the peptide and the m/z values

    """

    peptides_mz_value = {}

    if charge == 0: 
        raise ValueError (f"Charge is 0. Division by 0 is not possible")
    if proton_mass < 0: 
        raise ValueError (f"Proton mass is negative")
    
    for pep in peptide_mass_map: 
        if peptide_mass_map[pep] <= 0: 
            raise ValueError (f"Peptide mass is negative")
        mz_value = (peptide_mass_map[pep] + charge*proton_mass) / charge
        peptides_mz_value[pep] = mz_value
    
    return peptides_mz_value


def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42):
    """
    Plots a mass spectrum as a bar chart with m/z values on the x-axis and random intensities on the y-axis

    Parameters
    ----------
    mz_values : list
        list containing the m/z values for the spectrum
    random_count_range : tuple of int, optional
        Minimum and maximum values for randomly generated intensities with the defaultset to (0, 30000)
    seed : int, optional
        Random seed for reproducibility with the default set to 42

    Returns
    -------
    None : barchart
        bar chart with m/z values on the x-axis and random intensities on the y-axis
    """

    np.random.seed(seed)

    intensities = np.random.randint(random_count_range[0], random_count_range[1] + 1, size=len(mz_values)) # +1 so that the max is still included

    plt.figure(figsize=(10, 6))
    plt.bar(mz_values, intensities, width=0.5, color='lightgreen', edgecolor='black')
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.title("Simulated Mass Spectrum")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()


def fragment_peptide(peptide):
    """
    function to combine the b-ion and y-ion sequences for the peptide 

    Parameters
    ----------
    peptide: str   
        aminoacid sequence of the peptide

    Returns
    peptide_fragment : list
        list containing the combined b-ion ans y-ion series
    -------

    """
    b_ions = [peptide[:i] for i in range(1, len(peptide))] # cut from start to i; make sure whole peptide is not in it
    y_ions = [peptide[i:] for i in range(1, len(peptide))] # cut from i to end

    peptide_fragment = b_ions + y_ions
 
    return peptide_fragment