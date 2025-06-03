"""
module1.py

This module contains all the functions used to extract the coverage values from a dataframe
and to calculate the ratios.
Author: Gionmattia Carancini
Date: 07-03-2025
"""

# Import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def extract_single_region(df, study, coordinates, offset_up = 0, offset_down = 0):
    '''
    Code to extract a single region of interest for a single study, from the dataframe
    '''

    study_df = df[[study]].copy()
    # Question: Do I want to keep the offset here or do I want to let it be "user-defined"?
    region_df = study_df.iloc[coordinates["region_n"][0]:coordinates["region_n"][1]]

    return region_df


def compute_ratio(df, coordinates, N = 1000, random_state = None):
    '''
    Computes the ratio between the two selected regions. Requires slices of the three regions
    coordinates object must be a dictionary with following structure:
    coordinates = {x_slice: (a,b), y_slice: (a,b), n_slice: (a,b)}
    '''

    # Extract regions y and n
    region_y = df.iloc[coordinates["y_slice"][0]: coordinates["y_slice"][1]]
    region_n = df.iloc[coordinates["n_slice"][0]: coordinates["n_slice"][1]]
    # Extract and averages region x
    avg_x = df.iloc[coordinates["x_slice"][0]: coordinates["x_slice"][1]].mean()

    # averages regions y and n
    avg_y = region_y.mean()
    avg_n = region_n.mean()

    # calculates ratio translational activity (rta)
    rta = (avg_y - avg_n)/(avg_x - avg_n)

    # Convert to arrays, since it's faster
    y_array = region_y.to_numpy()
    n_array = region_n.to_numpy()

    combined = np.concatenate((y_array, n_array), axis=0)
    # get size y
    y_size = y_array.shape[0]

    # Prepare arrays to store permutation averages for each column
    perm_avg_y = np.empty((N, combined.shape[1]))
    perm_avg_n = np.empty((N, combined.shape[1]))

    rng = np.random.default_rng(random_state)

    for i in range(0, N):
        # Shuffle within each column
        # np.apply_along_axis applies a function to each column (axis=0)
        shuffled = np.apply_along_axis(rng.permutation, 0, combined)
        # Compute the means: first n_y rows for y, remaining for n
        perm_avg_y[i, :] = shuffled[:y_size, :].mean(axis=0)
        perm_avg_n[i, :] = shuffled[y_size:, :].mean(axis=0)

    fake_rta_ratios = (perm_avg_y - perm_avg_n) / (avg_x - perm_avg_n)

    # Compute p-values for each column (vectorised like here? or just dataframe?)
    # OLD version
    # p_values = np.mean(fake_ratio_dist >= rta, axis=0)
    p_values = np.mean(np.abs(fake_rta_ratios) >= np.abs(rta), axis = 0)

    # Convert to a DataFrame for clarity
    p_values_df = pd.DataFrame({
        'p_value': p_values
    }, index=df.columns)  # Assuming df has the original study columns

    print(p_values_df)

# Create a new dataframe with the results.
    result_df = pd.DataFrame({
        'avg_x': avg_x,
        'avg_y': avg_y,
        'avg_n': avg_n,
        'rta': rta
    })

    result_df["p_value"] = p_values

    return result_df


def cov_norm_df(raw_df):
    '''
    Computes the coverage as in (total raw counts on transcript) and divides each position by that number.
    '''

    col_sums = raw_df.sum()
    df = raw_df.div(col_sums, axis =1)

    return df


def zscore_norm_df(raw_df):
    '''
    Z-score normalisation.
    '''

    df =  (raw_df - raw_df.mean())/raw_df.std()

    return df


# A function to load files from a folder and set up pathways quickly?
### NEED: I need to add info on plot: UTRs, main CDS... frameshift...
def agg_score_plot(data, cds_start, cds_stop, site):
    # Aggregate z-scores by row into a new column
    data["RowAvg"] = data.mean(axis=1)

    # Normalize values for colormap
    norm = plt.Normalize(data["RowAvg"].min(), data["RowAvg"].max())
    colors = plt.cm.rainbow(norm(data["RowAvg"]))

    # Plot
    fig, ax = plt.subplots(figsize=(8, 1))
    ax.imshow([colors], aspect="auto", extent=[0, len(data), 0, 1])

    ax.set_xticks(np.arange(len(data)) + 0.5)
    ax.set_xticklabels(data.index)
    ax.axvline(x=cds_start, color='black', linestyle='--', linewidth=1.5)
    ax.axvline(x=cds_stop, color='black', linestyle='--', linewidth=1.5)
    ax.axvline(x=site, color='black', linestyle='--', linewidth=1.5)
    ax.set_yticks([])
    ax.set_title("1D Heat-Bar of Row Sums")

    plt.show()
    return







