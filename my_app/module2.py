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
def _clean_data(df, verbose=False):
    """
    Clean input dataframe by:
      1) Dropping columns that are all zeros.

    Args:
        df (pd.DataFrame): input data.
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


def _compute_xn_diffs(df, coordinates):
    """
    Lorem Ipsum

    """
    # slice out each region
    reg_x = df.iloc[coordinates["x_slice"][0]: coordinates["x_slice"][1]]
    reg_n = df.iloc[coordinates["n_slice"][0]: coordinates["n_slice"][1]]

    # Calculated as it follows so that the Pandas series keeps the rownames
    x_avg = reg_x.mean()
    n_avg = reg_n.mean()

    # Compute the absolute difference between x and n
    abs_xn_dff = abs(x_avg - n_avg)

    # Stores the values into a dataframe

    diffs_df = pd.DataFrame({
        "x avg": x_avg,
        "n avg": n_avg,
        "abs_xn_diff": abs_xn_dff}
    )

    return diffs_df, reg_x.to_numpy(), reg_n.to_numpy()


def _permute_values(base_array, noise_array, N, random_state=None):  # removed y_array, not needed.
    """
    Perform N permutations by shuffling combined baseline+noise arrays,
    recomputing fake difference each time.

    :param base_array: array with values from the "baseline" region of the permutation analysis
    :param noise_array: array with values from the "noise" region.
    :return perms_array: Returns an (N, n_series)-shaped array of "fake differences" for the permutation test.
    """

    # Combine the two arrays
    combined = np.concatenate((base_array, noise_array), axis=0)
    # Retrieve the original size of x
    base_size = base_array.shape[0]
    # Set the random state for permutations
    rng = np.random.default_rng(random_state)

    # Set the empty array for the results
    perms_array = np.empty((N, base_array.shape[1]))

    for i in range(N):
        # Shuffles the labels
        shuffled = np.apply_along_axis(rng.permutation, 0, combined)
        # Gets permutated sets averages
        perm_base = shuffled[:base_size].mean(axis=0)
        perm_noise = shuffled[base_size:].mean(axis=0)
        # Computes the diff between perm x and n. Stores it in the output array.
        perms_array[i, :] = abs(perm_base - perm_noise)

    return perms_array

# Computes the p-value and confidence boundary
def test_x_region(diffs_df, perms_array, ci=95):
    """
    One-tailed p-values and upper confidence bound for each series:
      p = (count(perms_array >= observed) + 1) / (N + 1)
      ci_upper = percentile(perms_array, ci)
    """

    N = perms_array.shape[0]
    # New code, we extract the observed difference between x and n
    observed = diffs_df["abs_xn_diff"]
    # Now we observe how many times the perms_array abs(x-n) were more extreme than the observed values.
    greater = (perms_array >= observed.values).sum(axis=0)
    pvals = (greater + 0.5) / (N + 1)

    #alpha = 100 - ci
    #lower = np.percentile(perms_array, alpha/2, axis = 0)
    #upper = np.percentile(perms_array, 100 - (alpha/2), axis=0)
    reject, pvals_bh, _, _ = multipletests(pvals, alpha= 1 - (ci/100), method='fdr_bh')

    return pd.DataFrame({
        "x avg": diffs_df["x avg"],
        "n avg": diffs_df["n avg"],
        #f'ci_lower_{ci}%': lower, # Not sure if I need it
        #f'ci_upper_{ci}%': upper, # Not sure if I need it
        'padj_bh': pvals_bh,
        'reject' : reject,
    }, index=observed.index)


def _compute_observed_ratios(valid_runs_df: pd.DataFrame, coordinates: dict[str, tuple[int, int]], good_runs: list
                             ) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract the three regions (x, y, and n)
    Compute observed ratio of translational-activity (RTA):
      (mean(region_y) - mean(region_n)) / (mean(region_x) - mean(region_n))

    Returns the ratios, and the two regions as numpy arrays.

    Args:
       valid_runs_df (pd.DataFrame): input data.A dataframe representing the coverage per nt position on a transcript. Each row is
        a position, each column is a different SRR (Run), among those that were prefiltered.
        coordinates: a dictionary, defining the boundaries of the regions "x_slice", "y_slice", and "n_slice". User-given

    Returns:
        ratios_df (pd.DataFrame): dataframe wtih averages for each slice, RTA.
        reg_y (np.ndarray): array of the values from y slice
        reg_n (np.ndarray): array of the values from n slice

    """

    # Keep only columns in "good_runs"
    df = valid_runs_df[good_runs].copy()

    # slice out each region
    reg_x = df.iloc[coordinates["x_slice"][0]: coordinates["x_slice"][1]]
    reg_y = df.iloc[coordinates["y_slice"][0]: coordinates["y_slice"][1]]
    reg_n = df.iloc[coordinates["n_slice"][0]: coordinates["n_slice"][1]]

    # Calculated as it follows so that the Pandas series keeps the rownames
    x_avg = reg_x.mean()
    y_avg = reg_y.mean()
    n_avg = reg_n.mean()

    # Compute the observed RTA
    RTA = (y_avg - n_avg) / (x_avg - n_avg)

    # Stores the values into a dataframe
    ratios_df = pd.DataFrame({
        "x avg": x_avg,
        "y avg": y_avg,
        "n avg": n_avg,
        "RTA": RTA,
    })

    # We keep the average x to calculate the fake ratios during the permutation test
    x_avg_array = x_avg.to_numpy()

    return ratios_df, reg_y.to_numpy(), reg_n.to_numpy(), x_avg_array


def _permute_ratios(y_array: np.ndarray,
                    n_array: np.ndarray,
                    x_avg_array: np.ndarray,
                    N: int,
                    random_state: bool = None
                    ) -> np.ndarray:  # removed y_array, not needed.
    """
    Perform N permutations by shuffling combined baseline+noise arrays,
    recomputing fake difference each time.

    :param base_array: array with values from the "baseline" region of the permutation analysis
    :param noise_array: array with values from the "noise" region.
    :return perms_array: Returns an (N, n_series)-shaped array of "fake differences" for the permutation test.
    """

    # Combine the two arrays
    combined = np.concatenate((y_array, n_array), axis=0)
    # Retrieve the original size of x
    y_size = y_array.shape[0]
    # Set the random state for permutations
    rng = np.random.default_rng(random_state)

    # Set the empty array for the results
    perms_array = np.empty((N, y_array.shape[1]))

    for i in range(N):
        # Shuffles the labels
        shuffled = np.apply_along_axis(rng.permutation, 0, combined)
        # Gets permutated sets averages
        perm_y = shuffled[:y_size].mean(axis=0)
        perm_n = shuffled[y_size:].mean(axis=0)
        # Computes the fake ratios for this permutation and stores them in the output array
        perms_array[i, :] = (perm_y - perm_n) / (x_avg_array - perm_n)

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
    observed = ratios_df["RTA"]
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
        'log2RTA': np.log2(ratios_df["RTA"]),
        'p-value': pvals,
        f'ci_lower_{ci}%': lower,
        f'ci_upper_{ci}%': upper,
        'padj_bh': pvals_bh,
        'significant' : reject,
    }, index=observed.index)


# Function linking together the various bits to be run all together
def run_analysis(df: pd.DataFrame,
                 coordinates: dict[str, tuple[int, int]],
                 ci: int,
                 N: int = 1000,
                 n_jobs: int = 3
                 ) -> pd.DataFrame:
    # Remove 0-count studies
    df = _clean_data(df)

    # Compute ratios and obtain arrays
    diffs_df, x_array, n_array = _compute_xn_diffs(df, coordinates)

    # Run in parallel
    perms_per_job = N // n_jobs
    results = Parallel(n_jobs=n_jobs)(
        delayed(_permute_values)(x_array, n_array, perms_per_job, random_state=i)
        for i in range(n_jobs)
    )

    # Stack the results into a single array (N, features)
    perms_array = np.vstack(results)

    # Compile dataframe of valid SRRs (they passed the perm test)
    valid_runs_df = test_x_region(diffs_df, perms_array, ci=ci)

    # Extract significant SRR from valid_runs_df
    good_runs = valid_runs_df.index[valid_runs_df["reject"] == True].tolist()

    # Function to open original dataframe, remove the rows not needed, calculate ratios with permutations...
    ratios_df, y_array, n_array, x_avg_array = _compute_observed_ratios(df, coordinates, good_runs)

    # Do the permutation test again, but only for the "valid" SRRs. In parallel.
    perms_per_job = N // n_jobs
    results = Parallel(n_jobs=n_jobs)(
        delayed(_permute_ratios)(y_array, n_array, x_avg_array, perms_per_job, random_state=i)
        for i in range(n_jobs)
    )
    # Stack the results into a single array (N, features)
    perms_array = np.vstack(results)

    # Attach values to the results_df, remove/rename values for x vs n...
    results_df = compute_pvalues_and_ci(ratios_df, perms_array, ci=95)

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

