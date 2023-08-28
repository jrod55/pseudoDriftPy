"""
Porting simulate_data.R functions
"""

import numpy as np
import pandas as pd

""" Functions for each batch mode """

""" Monotonic """
def t1(df, seed, m_eff):
    np.random.seed(seed)
    def calculate_area(row):
        if row['up_down']:
            return row['area'] * row['up']
        else:
            return row['area'] * row['down']

    # Equivalent to mutate in R
    df['up_down'] = np.random.choice([True, False], size=len(df))
    df['up'] = np.sort(np.absolute(np.random.uniform(low=1, high=m_eff, size=len(df))))
    df['down'] = np.sort(np.absolute(np.random.uniform(low=1, high=m_eff, size=len(df))))[::-1]
    df['area'] = df.apply(calculate_area, axis=1)
    
    # Equivalent to select in R
    df = df.drop(columns=['up_down', 'up', 'down'])
    
    return df

## Dummy data
data = {
    'batch': [1, 1, 2, 2, 3, 3],
    'area': [10, 15, 20, 25, 30, 35]
}
x = pd.DataFrame(data)

## Testing out the function with dummy data
seed = 42
m_eff = 10
xt1 = t1(x, seed, m_eff)

## Print results
print(x)
print(xt1)