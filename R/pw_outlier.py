import pandas as pd
import numpy as np
from scipy.spatial import distance
from scipy.stats import rankdata, norm, mad
import seaborn as sns
import matplotlib.pyplot as plt

def pw_outlier(df=None, n_cores=1, mad_threshold=3, pw_threshold=0.95, 
               peak_shrinkage=True, grouping_factor="batch", return_plot=False, 
               plot_name="pw_outlier_plot", samps_exclude="QC"):
    
    # Calculate most differing observations
    def get_most_diff(pairwise_diffs, current_geno, threshold_prop=0.34):
        geno_diffs = pd.DataFrame(pairwise_diffs, columns=current_geno['uid'].values)
        geno_diffs['sum_diffs'] = geno_diffs.sum(axis=1)
        top_diff_indices = geno_diffs['sum_diffs'].nlargest(int(len(geno_diffs) * threshold_prop)).index
        most_diff = current_geno.iloc[top_diff_indices].copy()
        most_diff['sum_diffs'] = geno_diffs['sum_diffs'].iloc[top_diff_indices].values
        return most_diff
    
    # STILL TO DO:
    # For function above, dictate threshold_prop in main function 
    # Check input df colnames
    # Identify most differing observations based on pairwise differences above the threshold
    # Include plotting capability to the function and the ability to save the plots to output
