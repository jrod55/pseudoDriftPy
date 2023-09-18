import pandas as pd
import numpy as np
from scipy.stats import rankdata

def pseudo_sdc(df=None,
               n_cores=1,
               train_batch=None,
               test_breaks=None, 
               test_window=None,
               test_index=None,
               criteria="MSE",
               qc_label=None,
               qc_multibatch=False,
               min_qc=5,
               quantile_increment=1,
               log_transform=True,
               mad_outlier=True,
               mad_threshold=3):
    
    nc = list(range(1, len(df['compound'].unique()) + 1))
    
    if len(nc) > 1:
        print(f"df contains more than one compound. Running for {max(nc)} compounds")
    
    comps_within = df['compound'].unique().tolist()

    # QCRSC helper function for Quality Control-Robust Spline Correction 
    def m_qcrsc(x, y, r):
        t_meta = [col for col in x.columns if col not in ["name", "compound", "area"]]
        tqc = x.pivot(index='compound', columns='name', values='area')
        
        # Placeholder for pmp::QCRSC function
        mc = None  # TODO: Need to implement a spline correction
        
        z = x.copy()
        z['area_corrected'] = mc[1]  # Placeholder
        
        m_qc = z[z['class'] == "QC"].copy()
        m_qc['rsd_robust'] = m_qc['area'].mad() / abs(m_qc['area'].median())
        m_qc['rsd_tqc_robust'] = m_qc['area_corrected'].mad() / abs(m_qc['area_corrected'].median())
        m_qc['rsd'] = m_qc['area'].std() / abs(m_qc['area'].mean())
        m_qc['rsd_tqc'] = m_qc['area_corrected'].std() / abs(m_qc['area_corrected'].mean())
        
        if r == "yes":
            return z
        else:
            return m_qc
    
    # Helper function
    def fun_within(x):
        sub_df = df[df['compound'] == x].copy()
        sub_df['class'] = np.where(sub_df['sample'].str.contains(qc_label), "QC", "Sample")
        sub_df['my_mad'] = sub_df.groupby(['batch', 'compound'])['area'].transform(lambda x: x.mad())
        sub_df['median_area'] = sub_df.groupby(['batch', 'compound'])['area'].transform('median')
        sub_df['high_thresh'] = sub_df['median_area'] + mad_threshold * sub_df['my_mad']
        sub_df['low_thresh'] = sub_df['median_area'] - mad_threshold * sub_df['my_mad']
        sub_df['area_tmp'] = np.where((sub_df['area'] >= sub_df['high_thresh']) | (sub_df['area'] <= sub_df['low_thresh']), np.nan, sub_df['area'])
        sub_df['area_imputed'] = sub_df.groupby(['batch', 'compound'])['area_tmp'].transform(lambda x: x.interpolate(method='linear', limit_area='inside'))
        m = sub_df[sub_df['batch'].isin(train_batch)].copy()
        m['my_mad'] = m.groupby(['batch', 'compound'])['area'].transform(lambda x: x.mad())
        m['median_area'] = m.groupby(['batch', 'compound'])['area'].transform('median')
        m['high_thresh'] = m['median_area'] + mad_threshold * m['my_mad']
        m['low_thresh'] = m['median_area'] - mad_threshold * m['my_mad']
        m['area_tmp'] = np.where((m['area'] >= m['high_thresh']) | (m['area'] <= m['low_thresh']), np.nan, m['area'])
        
        if mad_outlier:
            m['area_og'] = m['area']
            m['area'] = m['area_tmp']
            sub_df['area_og'] = sub_df['area']
            sub_df['area'] = sub_df['area_tmp']
        else:
            m['area_og'] = m['area']
            sub_df['area_og'] = sub_df['area']
        
        # Calculate vals_keep
        vals_keep = m[m['class'] == 'QC'].copy()
        vals_keep['pool_rank'] = vals_keep['area'].rank(pct=True)
        vals_keep = np.ceil(vals_keep['pool_rank'].agg(['min', 'max']) / quantile_increment) * quantile_increment
        
        # Calculate test_low, test_high, n_breaks, and k
        test_low = np.arange(0, vals_keep[0], quantile_increment)
        test_high = np.arange(vals_keep[1], 1, quantile_increment)
        n_breaks = test_breaks
        k = test_window
        
        # Placeholder for the dat_portion calculations per batch
        # TODO: Translate the dat_portion part of the function..what portion of the data and what partitioning of the batch gives best fit for the training data
        
        return sub_df  # Placeholder
    
    # TODO: Continue translating to apply QC-RSC to train set then test set
