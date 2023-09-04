import pandas as pd
import numpy as np

## Generate dummy data to test functions
def create_data():
    np.random.seed(seed)
    batches = ["batch_1", "batch_2", "batch_3"]
    sizes = [100, 200, 300]
    sample_names = []
    batches_column = []
    areas = []
    batch_indices = []
    
    for b, size in zip(batches, sizes):
        current_samples = []
        current_indices = []
        for i in range(size):
            if (i + 1) % 25 == 0:
                current_samples.append("QC")
            else:
                current_samples.append(f"sample_{i+1}")
            current_indices.append(i + 1)
        sample_names.extend(current_samples)
        batches_column.extend([b] * size)
        areas.extend(np.random.uniform(10, 50, size))
        batch_indices.extend(current_indices)

    df = pd.DataFrame({
        "sample": sample_names,
        "batch": batches_column,
        "batch_index": batch_indices,
        "area": areas
        
    })
    
    return df

## Input parameters to the functions and the dummy data
seed = 1234
m_eff = 10
b_eff_pos = 1.25
b_eff_neg = 0.75
x = create_data()

""" Functions for simulating each type of signal drift encountered in metabolomics experiments"""

## Montotonic increasing or decreasing in each batch
def t1(x):
    x = x.copy()
    np.random.seed(seed)
    
    # Helper function for manipulating by group
    def transform_group(group):
        n = len(group)
        up_down = np.random.choice([True, False])
        up = sorted(np.abs(np.random.uniform(1, m_eff, n)))
        down = sorted(np.abs(np.random.uniform(1, m_eff, n)), reverse=True)
        if up_down:
            group['area'] = group['area'] * up
        else:
            group['area'] = group['area'] * down
        return group[['sample', 'batch', 'area']]
    
    # Apply the helper function to each group and then concatenate the results
    return x.groupby('batch').apply(transform_group).reset_index(drop=True)

## Batch to batch
def t2(x):
    x = x.copy()
    np.random.seed(seed)
    
    def transform_group(group):
        blk = np.random.uniform(b_eff_neg, b_eff_pos)
        group['area'] = group['area'] * blk
        return group[['sample', 'batch', 'area']]
    
    return x.groupby('batch').apply(transform_group).reset_index(drop=True)

## Random
def t3(x):
    x = x.copy()
    np.random.seed(seed)
    
    def transform_group(group):
        n = len(group)
        group['area'] = group['area'] * np.abs(np.random.uniform(size=n))
        return group

    return x.groupby('batch').apply(transform_group).reset_index(drop=True)

def rleid(series):
    return (series != series.shift()).cumsum()

## Monotonic and batch-to-batch
def t4(x):
    x = x.copy()
    np.random.seed(seed)
    mm_eff = sorted([m_eff, 1])
    bb_eff = sorted([b_eff_neg, b_eff_pos])
    
    def transform_area(group):
        n = len(group)
        up_down = np.random.choice([True, False])
        up = sorted(np.abs(np.random.uniform(mm_eff[0], mm_eff[1], n)))
        down = sorted(np.abs(np.random.uniform(mm_eff[0], mm_eff[1], n)), reverse=True)
        group['area'] = np.where(up_down, group['area'] * up, group['area'] * down)
        return group
    
    def transform_blk(group):
        blk = np.random.uniform(bb_eff[0], bb_eff[1])
        group['area'] = group['area'] * blk
        return group

    qc_count = (x["sample"] == "QC").sum()
    x = x.assign(n_tile = x.groupby("batch")["batch_index"].transform(lambda y: pd.qcut(y, 2 * qc_count, labels=False, duplicates='drop')))
    x = x.groupby("batch").apply(transform_area).reset_index(drop=True)

    x = x.assign(up_down1 = x.groupby(["batch", "n_tile"])["batch_index"].transform(lambda _: np.random.choice([True, False])))
    x = x.assign(n_tile1 = x.groupby("batch")["up_down1"].transform(rleid))
    
    x = x.groupby(["batch", "n_tile1"]).apply(transform_blk).reset_index(drop=True)
    x_out = x.drop(columns=['up_down1', 'n_tile', 'n_tile1'])
    
    return x_out


## Test the functions
print(x.head())
xt1 = t1(x)
print(xt1.head())
xt2 = t2(x)
print(xt2.head())
xt3 = t3(x)
print(xt3.head())
xt4 = t4(x)
print(xt4.head())