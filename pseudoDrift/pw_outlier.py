import pandas as pd
import numpy as np
from scipy.spatial import distance
from scipy.stats import rankdata, norm
from scipy.stats import median_abs_deviation as mad
import seaborn as sns
import matplotlib.pyplot as plt

# Function to compute median absolute deviation
def compute_mad(arr):
    med = np.median(arr)
    return np.median(np.abs(arr - med))

## Main function
def pw_outlier(df=None, n_cores=1, mad_threshold=3, pw_threshold=0.95, 
               peak_shrinkage=True, grouping_factor="batch", samps_exclude="QC"):
    
    # Helper function to calculate most differing observations
    def get_most_diff(pairwise_diffs, current_geno, threshold_prop=0.05):
        num_samples = current_geno.shape[0]
        geno_diffs = pd.DataFrame(pairwise_diffs.reshape((num_samples, num_samples)), columns=current_geno['uid'].values, index=current_geno['uid'].values)
        geno_diffs['sum_diffs'] = geno_diffs.sum(axis=1)
        top_diff_indices = geno_diffs['sum_diffs'].nlargest(int(len(geno_diffs) * threshold_prop)).index
        most_diff = current_geno[current_geno['uid'].isin(top_diff_indices)].copy()
        most_diff['sum_diffs'] = geno_diffs['sum_diffs'].loc[top_diff_indices].values
        return most_diff
    
    # Check input df colnames
    c_names = df.columns
    c_names_need = ["name", "sample", "batch", "compound", "area", "rep", "rep_tech"]
    c_names_check = [col for col in c_names_need if col not in c_names]
    
    if len(c_names_check) >= 1:
        raise ValueError(f"Please check your input df to make sure it contains all necessary columns... {c_names_check} missing")
    
    m = df.copy()
    m[grouping_factor] = m[grouping_factor].astype('category')
    m['g'] = m[grouping_factor].cat.codes
    
    # Additional df calculations by compound
    m['geno_tmp'] = "TR" + m['rep_tech'].astype(str) + "_" + m['sample'].astype(str)
    m['tmp_rep'] = m.groupby(['geno_tmp', 'compound', 'rep_tech', 'g']).cumcount() + 1
    m['uid'] = m['geno_tmp'] + "_" + m['tmp_rep'].astype(str) + "_" + m['g'].astype(str) + "_" + m['compound']
    
    grouped = m.groupby(['g', 'compound'])
    m['my_rank'] = grouped['area'].transform(lambda x: rankdata(x))
    m['my_quant'] = norm.ppf(m['my_rank'] / (grouped['area'].transform('size') + 1))
    m['my_mad'] = grouped['area'].transform(lambda x: compute_mad(x))
    m['my_thresh'] = mad_threshold * m['my_mad'] + grouped['area'].transform('median')
    
    m['area_tmp'] = np.where(m['area'] > m['my_thresh'], np.nan, m['area'])
    m['area_shrink'] = np.where(m['area'] > m['my_thresh'], grouped['area_tmp'].transform('max') + m['my_quant'], m['area'])
    
    if peak_shrinkage:
        m['area_og'] = m['area']
        m['area'] = m['area_shrink']
    else:
        m['area_og'] = m['area']
        m['area'] = m['area_tmp']
    
    m['uid_compound'] = m['uid'] + "_" + m['compound']
    mog = m.copy()
    m = m[~m['sample'].isin([samps_exclude])]
    
    # Compounds in data
    integrated_compounds = m['compound'].unique()
    samps_rm = []
    samps_plot = []
    
    # Loop through compounds
    for compound in integrated_compounds:
        my_df = m[m['compound'] == compound]
        my_batches = my_df['g'].unique()
        diff_vect = []
        most_diff = []
        plot_dat = []
        
        for batch in my_batches:
            current_batch = my_df[my_df['g'] == batch]
            my_genos = current_batch['geno_tmp'].unique()
            diff_vect1 = []
            most_diff1 = []
            
            for geno in my_genos:
                current_geno = current_batch[current_batch['geno_tmp'] == geno]
                current_geno_name = geno
                
                if len(current_geno) < 3:
                    most_diff1.append(None)
                    diff_vect1.append(None)
                else:
                    pairwise_diffs = distance.squareform(distance.pdist(current_geno['area'].values.reshape(-1, 1), 'cityblock'))
                    most_diff1.append(get_most_diff(pairwise_diffs, current_geno))
                    diff_vect1.append(pd.DataFrame({
                        'value': pairwise_diffs.flatten(),
                        'compound': compound,
                        'g': batch
                    }))
            
            if len(diff_vect1) == 0:
                diff_vect.append(None)
                plot_dat.append(None)
                most_diff.append(None)
            else:
                diff_vect.append(pd.concat(diff_vect1))
                most_diff.append(pd.concat(most_diff1))

        samps_rm.append(pd.concat(most_diff))
        samps_plot.append(pd.concat(diff_vect))
    
    # Calculate the threshold for each batch and compound in each batch
    thresholds = []
    for df in samps_plot:
        if df is not None:
            thresh = df.groupby(['compound', 'g'])['value'].quantile(pw_threshold).reset_index()
            thresholds.extend(thresh.to_dict(orient='records'))
    
    # Filter out outliers
    outliers = pd.concat(samps_rm)
    outliers = outliers.drop_duplicates(subset=['uid', 'compound'])
    outliers['uid_compound'] = outliers['uid'] + "_" + outliers['compound']
    m1 = mog[~mog['uid_compound'].isin(outliers['uid_compound'])]
    m1['area'] = m1['area_og']
    m1 = m1[c_names]
    
    return {
        'df_cleaned': m1,
        'df_rm': outliers[c_names],
        'thresholds': thresholds
    }


# Generate larger sample dataset to test function
n_samples = 1000
batches = ['batch_' + str(i) for i in range(1, 6)]
compounds = ['compound_' + str(i) for i in range(1, 4)]
reps = [1, 2, 3]
rep_techs = ['tech_' + str(i) for i in range(1, 4)]

df_sample = pd.DataFrame({
    'name': ['name_' + str(i) for i in range(n_samples) for _ in compounds for _ in reps for _ in rep_techs],
    'sample': ['sample_' + str(i) for i in range(n_samples) for _ in compounds for _ in reps for _ in rep_techs],
    'batch': np.random.choice(batches, n_samples*len(compounds)*len(reps)*len(rep_techs)),
    'compound': np.tile(np.repeat(compounds, len(reps)*len(rep_techs)), n_samples),
    'rep': np.tile(np.repeat(reps, len(rep_techs)), n_samples*len(compounds)),
    'rep_tech': np.tile(rep_techs, n_samples*len(compounds)*len(reps)),
})

# Simulate area values for each compound with some outliers
compound_means = {'compound_1': 100, 'compound_2': 200, 'compound_3': 50}
std_dev = 10

for compound, mean_value in compound_means.items():
    mask = df_sample['compound'] == compound
    df_sample.loc[mask, 'area'] = np.random.normal(mean_value, std_dev, mask.sum())

    # Introducing some outliers
    outlier_indices = np.random.choice(df_sample[mask].index, size=int(0.05 * mask.sum()), replace=False)
    df_sample.loc[outlier_indices, 'area'] = mean_value + np.random.choice([50, -50], len(outlier_indices))

df_sample.head()

# Test the function
my_results = pw_outlier(df_sample)

## Notes: There is a slight bug that needs to be fixed in plotting function below

""" # Function for plotting
def generate_and_display_plots(df, thresholds):
    integrated_compounds = df['compound'].unique().tolist()

    for compound in integrated_compounds:
        compound_diffs = df[df['compound'] == compound]['area'].tolist()
        
        # Plot pairwise differences for the compound
        plt.figure(figsize=(10, 5))
        sns.kdeplot(compound_diffs, shade=True, color='lightblue')
        for batch_info in thresholds:
            if batch_info['compound'] == compound:
                plt.axvline(x=batch_info['threshold'], color='red', linestyle='--')
        plt.title(f"Compound: {compound}")
        plt.xlabel("Replication pairwise differences")
        plt.ylabel("Density")
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.tight_layout()
        plt.show()

## Test the plotting function
generate_and_display_plots(my_results['df_cleaned'], my_results['thresholds']) """
