# proteomics course package
- repository for the python advanced practical course

This is a proteomics simulation performing protein digestion, liquid chromatography and mass spectroscopy for protein with an enzyme. To start it needs a dictionary with the Proteins and their sequences, as well as the digestion sequence of the enzyme (if a different enzyme than the implemented enzymes is used).

The package enables the input of FASTA proteins, the digestion (producing the expected peptides and calculating the sequence coverage), liquid chromatography (predicting the retention times, plotting and filtering them) and the MS simulation (calculating the peptide mol masses and mz values, plotting the mz values and creating the fragments of the peptides).

The contained functions are the following: 

Fasta: 
- read_fasta (reads the FASTA file)
Protein digestion: 
- digest_protein_sequence (digests the protein with an enzyme and outputs the expected peptides)
- digest_protein_collection (digests multiple proteins with an enzyme and outputs the expected peptides)
- compute_sequence_coverage (computes the sequence coverage of the protein)

Liquid chromatography
- predict_lc_retention_times (predicts the retention times of the computed peptides)
- plot_retention_time (plots the retention times)
- select_retention_time_window (filters the peptides for peptides with retention times in a given window)

Mass spectra: 
- calculate_mol_mass (calculates the mol mass of one aminoacid sequence)
- calculate_mol_mass_collection (calculates the mol mass of multiple peptides)
- calculate_mz_collection (calculates the mz values for the peptides)
- plot_spectrum (plots the mz values of the peptides)
- fragment_peptide (gives the fragments for a given peptide sequence)

The workflow is demonstrated in the notebook ms_experiment_final.ipynb