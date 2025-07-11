"""
module2.py

This module contains all the functions used in RiboMARTA, as per the refactoring done in April 2025
Author: Gionmattia Carancini
Date: 23-04-2025
"""
import pandas
import pandas as pd
import numpy as np
from narwhals import Boolean
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
import matplotlib.pyplot as plt


# needed to remove samples where we do not have info
def _clean_data(df, coordinates, verbose=False):
    """
    Clean input dataframe by:
      1) Dropping columns that are all zeros.

    Args:
        df (pd.DataFrame): input data.
        coordinates (dict): must contain "x_slice" as (start, end) row indices.
        verbose (bool): if True, prints filtering stats.

    Returns:
        pd.DataFrame: cleaned dataframe.
    """
    # Drop all-zero columns
    original_shape = df.shape[1]
    df_clean = df.loc[:, (df != 0).any(axis=0)]
    zero_filtered_shape = df_clean.shape[1]

    # Originally, we used to remove all columns were x was 0...now we are doing the permutation test

    if verbose:
        print(f"Numbers of studies originally: {original_shape}")
        print(f"After removing studies with zero coverage: {zero_filtered_shape}")

    return df_clean


# Computers the ratios
def _compute_observed_ratios(df, coordinates):
    """
    Extract the three regions (x, y, and n)
    Compute observed ratio of translational-activity (RTA):
      (mean(region_y) - mean(region_n)) / (mean(region_x) - mean(region_n))

    Returns the ratios, and the three regions as numpy arrays.

    Args:
        df (pd.DataFrame): input data.
        coordinates (dict): must contain "x_slice" as (start, end) row indices.

    Returns:
        ratios_df (pd.DataFrame): dataframe wtih averages for each slice, RTA, and the absolute difference between average x and n.
        reg_x (np.ndarray): array of the values from x slice
        reg_n (np.ndarray): array of the values from n slice

    """
    # slice out each region
    reg_x = df.iloc[coordinates["x_slice"][0]: coordinates["x_slice"][1]]
    reg_y = df.iloc[coordinates["y_slice"][0]: coordinates["y_slice"][1]]
    reg_n = df.iloc[coordinates["n_slice"][0]: coordinates["n_slice"][1]]

    # Calculated as it follows so that the Pandas series keeps the rownames
    x_avg = reg_x.mean()
    y_avg = reg_y.mean()
    n_avg = reg_n.mean()

    RTA = (y_avg - n_avg) / (x_avg - n_avg)

    # Compute the absolute difference between x and n
    abs_xn_dff = abs(x_avg - n_avg)

    # Stores the values into a dataframe

    ratios_df = pd.DataFrame({
        "x avg": x_avg,
        "y avg": y_avg,
        "n avg": n_avg,
        "RTA": RTA,
        "abs_xn_diff": abs_xn_dff}
    )

    # y_avg_array = y_avg.to_numpy()

    return ratios_df, reg_x.to_numpy(), reg_n.to_numpy()


# Computes the permutations array
def _permute_ratios(x_array, n_array, N, random_state=None):  # removed y_array, not needed.
    """
    Perform N permutations by shuffling combined y+n,
    recomputing fake TA ratios each time.
    Returns an (N, n_series)-shaped array of fake ratios.
    """

    # Combine the two arrays
    combined = np.concatenate((x_array, n_array), axis=0)
    # Retrieve the original size of x
    x_size = x_array.shape[0]
    # Set the random state for permutations
    rng = np.random.default_rng(random_state)

    # Set the empty array for the results
    perms_array = np.empty((N, x_array.shape[1]))

    for i in range(N):
        # Shuffles the labels
        shuffled = np.apply_along_axis(rng.permutation, 0, combined)
        # Gets permutated sets averages
        perm_x = shuffled[:x_size].mean(axis=0)
        perm_n = shuffled[x_size:].mean(axis=0)
        # Computes the diff between perm x and n. Stores it in the output array.
        perms_array[i, :] = abs(perm_x - perm_n)

    return perms_array


# Computes the p-value and confidence boundary
def compute_pvalues_and_ci(ratios_df, permuted, ci=95):
    """
    One-tailed p-values and upper confidence bound for each series:
      p = (count(permuted >= observed) + 1) / (N + 1)
      ci_upper = percentile(permuted, ci)
    """

    N = permuted.shape[0]
    # New code, we extract the observed difference between x and n
    observed = ratios_df["abs_xn_diff"]
    # Now we observe how many times the permuted abs(x-n) were more extreme than the observed values.
    greater = (permuted >= observed.values).sum(axis=0)
    pvals = (greater + 0.5) / (N + 1)

    alpha = 100 - ci
    lower = np.percentile(permuted, alpha/2, axis = 0)
    upper = np.percentile(permuted, 100 - (alpha/2), axis=0)
    reject, pvals_bh, _, _ = multipletests(pvals, alpha= 1 - (ci/100), method='fdr_bh')

    return pd.DataFrame({
        "x avg": ratios_df["x avg"],
        "y avg": ratios_df["y avg"],
        "n avg": ratios_df["n avg"],
        'RTA': ratios_df["RTA"],
        'p-value': pvals,
        f'ci_lower_{ci}%': lower,
        f'ci_upper_{ci}%': upper,
        'padj_bh': pvals_bh,
        'reject' : reject,
    }, index=observed.index)


# Function linking together the various bits to be run all together
def run_analysis(df: pd.DataFrame,
                 coordinates: dict[str, tuple[int, int]],
                 ci: int,
                 N: int=1000,
                 n_jobs: int=3
                 ) -> pd.DataFrame:
    # NEED to FORMAT (remove columns) the dataframe

    # Remove 0-count studies
    df = _clean_data(df, coordinates)
    # Compute ratios and obtain arrays
    ratios_df, x_array, n_array = _compute_observed_ratios(df, coordinates)
    # Run in parallel
    perms_per_job = N // n_jobs
    results = Parallel(n_jobs=n_jobs)(
        delayed(_permute_ratios)(x_array, n_array, perms_per_job, random_state=i)
        for i in range(n_jobs)
    )
    # Stack the results into a single array (N, features)
    perms_array = np.vstack(results)
    # Compile results dataframe
    results_df = compute_pvalues_and_ci(ratios_df, perms_array, ci= ci)

    return results_df


# Function to integrate metadata with results
def integrate_metadata(results: pd.DataFrame,
                       metadata: pd.DataFrame
                       ) -> pd.DataFrame:
    """
    Merges RiboMARTA results with the available metadata, using "Run" ID as index.

    Args:
        results (pd.DataFrame): dataframe generated by run_analysis().
        metadata (pd.DataFrame): dataframe containing metadata, following the Riboseq.Org Data Format from Oct2024.

    Returns:
        pd.DataFrame: merged dataframe, with results and metadata for each SRR (Run).
    """
    # Creates a slim version of the metadata
    slim_meta = metadata[["Run", "TISSUE", "CELL_LINE", "Strain", "Genotype", "Growth_Condition", "Info", "SampleName", "sample_source", "BioSample"]].copy()
    # assigns the column "Run" as the index
    slim_meta.set_index("Run", inplace=True)
    # merge the two dataframes on index, keep only indexes in results
    merged_df =  pd.merge(slim_meta, results, left_index=True, right_index=True, how="right")

    return merged_df

