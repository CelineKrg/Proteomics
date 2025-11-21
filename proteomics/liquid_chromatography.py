from pyteomics import achrom
import matplotlib.pyplot as plt

def predict_lc_retention_times(peptides):
    """
    function to predict the retention time of the peptides

    Parameters
    peptides : list
        list containing all the peptides obtained after the digestion 
    ----------

    Returns
    retention_times : dict
        dictionary of the predicted retention times for the peptides
    -------
    """
    retention_times = {} # create empty dict

    for pep in peptides: 
        relat_RT = achrom.calculate_RT(pep, achrom.RCs_guo_ph7_0)
        retention_times[pep] = round(float(relat_RT), 2) # float to convert np.float into python float 

    return retention_times


def plot_retention_time(retention_times, resolution=30):
    """
    Plots a histogram of peptide retention times

    Parameters
    ----------
    retention_times : dict
        dictionary of retention times with peptide sequences as keys and retention times as values
    resolution : int, optional
        Number of bins in the histogram with the default set to 30

    Returns
    -------
    None : histogram
        histogram of the retention times
    """

    data = list(retention_times.values())

    plt.figure(figsize=(8, 5))
    plt.hist(data, bins=resolution, color='lightgreen', edgecolor='black')
    plt.xlabel("Relative retention time")
    plt.ylabel("Frequency")
    plt.title("Histogram of peptide retention times")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.show()


def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    function to filter the retention time of peptides within a given range

    Parameters
    ----------
    peptide_rt_map : dict
        dictionary of retention times with peptide sequences as keys and retention times as values

    lower_ret_time : int
        minimal retention time for filtering in minutes

    upper_ret_time : int
        maximal retention time for filtering in minutes
    
    Returns
    -------
    filtered_peptide_rt_map : dict
        dictionary containing peptides and their retention times within the given range

    """
    filtered_peptide_rt_map = {}

    for pep, rt in peptide_rt_map.items(): 
        if lower_ret_time <= rt <= upper_ret_time:
            filtered_peptide_rt_map[pep] = rt
    
    return filtered_peptide_rt_map