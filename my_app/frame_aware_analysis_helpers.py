"""
frame_aware_analysis_helpers.py

This module contains all the functions used in RiboMARTA, but readapted to work for a frame-aware analysis.
Only the functions readapted have been rewritten (and they have a "fa" in the name, as in "frame aware").
Other functions have been imported from module2.
Also, I added a couple of functions to prepare the data for plotting.

Author: Gionmattia Carancini
Date: 25-08-2025
"""
import pandas
import pandas as pd
import numpy as np
from narwhals import Boolean
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
from module2 import _clean_data
from module2 import _permute_values
from module2 import _permute_ratios
from module2 import compute_pvalues_and_ci

def _fa_compute_xn_diffs(df: pd.DataFrame, coordinates: dict, baseline_frame: int):
    """
    Extracts coverage info from df, using coordinates (they define baseline and noise region).
    Works only with the

    """
    # slice out each region
    slice_x = df.iloc[coordinates["x_slice"][0]: coordinates["x_slice"][1]]
    slice_n = df.iloc[coordinates["n_slice"][0]: coordinates["n_slice"][1]]

    # Subselect only the frame of interest, remove "Frame"
    reg_x = slice_x[slice_x["Frame"] == baseline_frame].copy() # THIS NEEDS TESTING
    reg_x.drop("Frame", axis=1, inplace=True)
    reg_n = slice_n.copy() # Playing conservative, we do not require a frame for noise region
    reg_n.drop("Frame", axis=1, inplace=True)

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

def _fa_compute_observed_ratios(clean_df: pd.DataFrame, coordinates: dict[str, tuple[int, int]], good_runs: list,
                             baseline_frame: int, target_frame: int) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract the three regions (x, y, and n)
    Compute observed ratio of translational-activity (RTA):
      (mean(region_y) - mean(region_n)) / (mean(region_x) - mean(region_n))

    Returns the ratios, and the two regions as numpy arrays.

    Args:
       clean_df (pd.DataFrame): input data.A dataframe representing the coverage per nt position on a transcript. Each row is
        a position, each column is a different SRR (Run), among those that were prefiltered.
        coordinates: a dictionary, defining the boundaries of the regions "x_slice", "y_slice", and "n_slice". User-given

    Returns:
        ratios_df (pd.DataFrame): dataframe wtih averages for each slice, RTA.
        reg_y (np.ndarray): array of the values from y slice
        reg_n (np.ndarray): array of the values from n slice

    """

    # Keep only columns in "good_runs"
    df = clean_df[good_runs].copy()
    # Add again frame info
    df["Frame"] = clean_df["Frame"]

    # slice out each region
    slice_x = df.iloc[coordinates["x_slice"][0]: coordinates["x_slice"][1]]
    slice_y = df.iloc[coordinates["y_slice"][0]: coordinates["y_slice"][1]]
    slice_n = df.iloc[coordinates["n_slice"][0]: coordinates["n_slice"][1]]

    # select only things in a specific frame
    reg_x = slice_x[slice_x["Frame"] == baseline_frame].copy()
    reg_x.drop("Frame", axis=1, inplace=True)
    reg_y = slice_y[slice_y["Frame"] == target_frame].copy()
    reg_y.drop("Frame", axis=1, inplace=True)
    reg_n = slice_n.copy()
    reg_n.drop("Frame", axis=1, inplace=True)

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


# Function linking together the various bits to be run all together
def run_fa_analysis(df: pd.DataFrame,
                 coordinates: dict[str, tuple[int, int]],
                 baseline_frame: int,
                 target_frames: list,
                 ci: int,
                 N: int = 1000,
                 n_jobs: int = 3
                 ) -> pd.DataFrame:

    # Remove 0-count studies
    df = _clean_data(df)

    # Compute ratios and obtain arrays
    diffs_df, x_array, n_array_1 = _fa_compute_xn_diffs(df, coordinates, baseline_frame)

    # Run in parallel
    perms_per_job = N // n_jobs
    results = Parallel(n_jobs=n_jobs)(
        delayed(_permute_values)(x_array, n_array_1, perms_per_job, random_state=i)
        for i in range(n_jobs)
    )

    # Stack the results into a single array (N, features)
    perms_array = np.vstack(results)

    # Compile dataframe of valid SRRs (they passed the perm test)
    valid_runs_df = test_x_region(diffs_df, perms_array, ci=ci)

    # Extract significant SRR from valid_runs_df
    good_runs = valid_runs_df.index[valid_runs_df["reject"] == True].tolist()

    # Set up dictionary of results
    results_dict = {}

    # Function to open original dataframe, remove the rows not needed, calculate ratios with permutations...
    for f in target_frames:
        ratios_df, y_array, n_array_2, x_avg_array = _fa_compute_observed_ratios(df,
                                                                            coordinates,
                                                                            good_runs,
                                                                            baseline_frame,
                                                                            f)

        # Do the permutation test again, but only for the "valid" SRRs. In parallel.
        perms_per_job = N // n_jobs
        results = Parallel(n_jobs=n_jobs)(
            delayed(_permute_ratios)(y_array, n_array_2, x_avg_array, perms_per_job, random_state=i)
            for i in range(n_jobs)
        )
        # Stack the results into a single array (N, features)
        perms_array = np.vstack(results)

        # Attach values to the results_df, remove/rename values for x vs n...
        results_df = compute_pvalues_and_ci(ratios_df, perms_array, ci=95)

        # Store df in results_dict
        results_dict[f"Frame_{f}"] = results_df

    return results_dict


def format_frame_aware_res(baseline_frame: int, results_dict: dict) -> pd.DataFrame:
    """
    Formats the results dictionary as a dataframe, for easiness of plotting
    :param baseline_frame: the frame for the baseline
    :param results_dict:  dictionary with the results of run_analysis
    :return:
    """

    results_df = pd.DataFrame()
    # Extract the series to compile into a single dataframe
    keys = list(results_dict.keys())

    srr_series = results_dict[keys[0]].index
    x_avg_series = results_dict[keys[0]]["x avg"]
    n_avg_series = results_dict[keys[0]]["n avg"]

    # Create a DataFrame from the Series directly
    results_df = pd.DataFrame({
        "SRR": srr_series,
        f"x avg (f{baseline_frame})": x_avg_series,
        "n avg": n_avg_series,
    })

    for k in keys:
        y_colname = f"y avg ({k})"
        RTA_colname = f"RTA ({k})"
        padj_colname = f"padj ({k})"
        results_df[y_colname] = results_dict[k]["y avg"]
        results_df[RTA_colname] = results_dict[k]["RTA"]
        results_df[padj_colname] = results_dict[k]["padj_bh"]

    return results_df

def setup_df_for_plot(results_df):
    """
    Formats the dataframe for the plotting function...

    :param results_df:
    :return: long_df
    """
    long_df = pd.melt(results_df,
                      id_vars=['SRR'],
                      value_vars=['RTA (Frame_1)', 'RTA (Frame_2)', 'RTA (Frame_3)'],
                      var_name='Frame',
                      value_name='RTA')
    long_df['Frame'] = long_df['Frame'].str.replace('RTA ', '')

    padj_df = pd.melt(results_df,
                      id_vars=['SRR'],
                      value_vars=['padj (Frame_1)', 'padj (Frame_2)', 'padj (Frame_3)'],
                      var_name='Frame',
                      value_name='padj')

    padj_df['Frame'] = padj_df['Frame'].str.replace('padj ', '')

    long_df = long_df.merge(padj_df[['SRR', 'Frame', 'padj']], on=['SRR', 'Frame'])

    long_df['Significance'] = long_df['padj'].apply(lambda x: 'Significant' if x < 0.05 else 'Not Significant')

    long_df['Frame'] = long_df['Frame'].str.replace('(', '')
    long_df['Frame'] = long_df['Frame'].str.replace(')', '')
    long_df['Frame'] = long_df['Frame'].str.replace('_', ' ')
    long_df['symbol_size'] = long_df['Significance'].apply(lambda x: 0.1 if x == 'Not Significant' else 0.5)

    long_df['hover_text'] = ('<b>SRR</b>: ' + long_df['SRR'].astype(str) + '<br>' +
                             '<b>RTA</b>: ' + long_df['RTA'].round(2).astype(str) + '<br>' +
                             '<b>Frame</b>: ' + long_df['Frame'].astype(str) + '<br>' +
                             '<b>padj</b>: ' + long_df['padj'].round(4).astype(str))

    return long_df


def RTA_by_frame_plotter(long_df):
    """
    Plotting Function.

    :param long_df:
    :return:
    """
    # Define custom colors and markers for the plot
    color_map = {'Frame 1': "#48cf4b", 'Frame 2': "#68ace3", 'Frame 3': "#e66343"}

    symbol_map = {'Significant': 'diamond', 'Not Significant': 'circle'}

    # Create the interactive plot using Plotly Express
    fig = px.scatter(long_df,
                     x='SRR',
                     y='RTA',
                     color='Frame',
                     symbol='Significance',
                     size='symbol_size',
                     color_discrete_map=color_map,
                     symbol_map=symbol_map,
                     hover_data=['padj', 'Frame'],
                     title='Interactive RTA Values per SRR by Frame and Significance')

    fig.update_traces(hovertemplate='<b>SRR</b>: %{x}<br>'
                                    '<b>RTA</b>: %{y:.2f}<br>'
                                    '<b>Frame</b>: %{customdata[1]}<br>'
                                    '<b>padj</b>: %{customdata[0]:.4f}<br>'
                                    '<extra></extra>')

    # Add vertical lines connecting the points
    for srr in long_df['SRR'].unique():
        fig.add_vrect(x0=srr, x1=srr,
                      y0=long_df[long_df['SRR'] == srr]['RTA'].min(),
                      y1=long_df[long_df['SRR'] == srr]['RTA'].max(),
                      line_width=2,
                      line_dash='solid',
                      line_color='black',
                      layer="below")

    # Update the layout for a cleaner look
    fig.update_layout(
        template="plotly_white",
        xaxis_title="SRR",
        yaxis_title="RTA"
    )

    return fig