"""
module2.py

This module contains all the functions used in RiboMARTA, as per the refactoring done in April 2025
Author: Gionmattia Carancini
Date: 23-04-2025
"""
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import numpy as np


# needed to remove samples where we do not have info
def _clean_data(df, coordinates, verbose=False):
    """
    Clean input dataframe by:
      1) Dropping columns that are all zeros.
      2) Dropping columns whose mean over the x-slice is zero.

    Args:
        df (pd.DataFrame): input data.
        coordinates (dict): must contain "x_slice" as (start, end) row indices.
        verbose (bool): if True, prints filtering stats.

    Returns:
        pd.DataFrame: cleaned dataframe.
    """
    # 1. Convert Nan to zero values
    df = df.fillna(0)

    # 2. drop all-zero columns
    original_shape = df.shape[1]
    df_filtered = df.loc[:, (df != 0).any(axis=0)]
    zero_filtered_shape = df_filtered.shape[1]

    # 3. drop columns with zero mean in x-slice
    x0, x1 = coordinates["x_slice"]
    avg_x = df_filtered.iloc[x0:x1].mean()
    zero_avg_cols = avg_x[avg_x == 0].index.tolist()
    df_clean = df_filtered.drop(columns=zero_avg_cols)
    zero_x_filtered_shape = df_clean.shape[1]

    if verbose:
        print(f"Numbers of studies originally: {original_shape}")
        print(f"After removing studies with zero coverage: {zero_filtered_shape}")
        print(f"After removing studies with no coverage in x region: {zero_x_filtered_shape}")

    return df_clean

# Computers the ratios
def _compute_observed_ratios(df, coordinates):
    """
    Extract the three regions (x, y, and n)
    Compute observed translational-activity ratios:
      (mean(region_y) - mean(region_n)) / (mean(region_x) - mean(region_n))

    Returns the ratios, and the three regioms as numpy arrays.
    """
    # slice out each region
    reg_x = df.iloc[coordinates["x_slice"][0]: coordinates["x_slice"][1]]
    reg_y = df.iloc[coordinates["y_slice"][0]: coordinates["y_slice"][1]]
    reg_n = df.iloc[coordinates["n_slice"][0]: coordinates["n_slice"][1]]

    # Calculated as it follows so that the Pandas series keeps the rownames
    x_avg = reg_x.mean()
    y_avg = reg_y.mean()
    n_avg = reg_n.mean()

    ratios = (y_avg - n_avg) / (x_avg - n_avg)

    # Stores the values intoa dataframe

    ratios_df = pd.DataFrame({
        "x avg": x_avg,
        "y avg": y_avg,
        "n avg": n_avg,
        "ratios": ratios}
    )

    x_avg_array = x_avg.to_numpy()

    return ratios_df, reg_x.to_numpy(), reg_y.to_numpy(), reg_n.to_numpy(), x_avg_array

# Computes the permutations array

def _permute_ratios(y_array, n_array, x_avg_array, N, random_state=None):
    """
    Perform N permutations by shuffling combined y+n,
    recomputing fake TA ratios each time.
    Returns an (N, n_series)-shaped array of fake ratios.
    """

    # Combine the two arrays
    combined = np.concatenate((y_array, n_array), axis=0)
    # Retrieve the original size of y
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
        perms_array[i,:] = (perm_y - perm_n)/ (x_avg_array - perm_n)

    return perms_array

# Computes the p-value and confidence boundary
def compute_pvalues_and_ci(ratios_df, permuted, ci=95):
    """
    One-tailed p-values and upper confidence bound for each series:
      p = (count(permuted >= observed) + 1) / (N + 1)
      ci_upper = percentile(permuted, ci)
    """

    N = permuted.shape[0]
    observed = ratios_df["ratios"]
    greater = (permuted >= observed.values).sum(axis=0)
    pvals = (greater + 0.5) / (N + 1)
    upper = np.percentile(permuted, ci, axis=0)
    reject, pvals_bh, _, _ = multipletests(pvals, alpha= 1 - (ci/100), method='fdr_bh')

    return pd.DataFrame({
        "x avg": ratios_df["x avg"],
        "y avg": ratios_df["y avg"],
        "n avg": ratios_df["n avg"],
        'observed_rta': observed.values,
        'p-value': pvals,
        f'ci_upper_{ci}%': upper,
        'padj_bh': pvals_bh,
        'significant' : reject,
    }, index=observed.index)

# Function linking together the various bits to be run all together
def run_analysis(df, coordinates, ci, N = 1000, n_jobs = 3):
    # NEED to FORMAT (remove columns) the dataframe

    # Remove 0-count studies
    df = _clean_data(df, coordinates)
    # Compute ratios and obtain arrays
    ratios_df, x_array, y_array, n_array, x_avg_array = _compute_observed_ratios(df, coordinates)
    # Run in parallel
    perms_per_job = N // n_jobs
    results = Parallel(n_jobs=n_jobs)(
        delayed(_permute_ratios)(y_array, n_array, x_avg_array, perms_per_job, random_state=i)
        for i in range(n_jobs)
    )
    # Stack the results into a single array (N, features)
    perms_array = np.vstack(results)
    # Compile results dataframe
    results_df = compute_pvalues_and_ci(ratios_df, perms_array, ci=95)

    return results_df

# Function to make the boxplot of RTAs
def plot_observed_rta_interactive(results_df,
                                  box_fill="gray",
                                  box_line_color="black",
                                  box_width=0.2,
                                  jitter=0.2):
    """
    Interactive box+points plot of “observed_rta” (where “significant” == True),
    with hover labels, custom box color, and controlled box width.
    """
    filtered = results_df.loc[results_df["significant"] == True, ["observed_rta"]].dropna()
    if filtered.empty:
        raise Warning("No significant RTAs were present among the results")

    y_vals = filtered["observed_rta"].values
    labels = filtered.index.astype(str).tolist()

    x_vals = [1] * len(y_vals)

    box = go.Box(
        x=x_vals,
        y=y_vals,
        boxpoints="all",
        jitter=jitter,
        pointpos=0,
        width=box_width,
        fillcolor= box_fill,
        line=dict(color=box_line_color),
        marker=dict(color="darkblue", opacity=0.6),
        text=labels,
        hoverinfo="text+y",
        name="",
        showlegend=False
    )

    layout = go.Layout(
        title=dict(
            text="Boxplot of Observed RTAs across studies",
            font=dict(color="black")
        ),
        xaxis=dict(
            showticklabels=False,
            range=[0.5, 1.5],
            showgrid = False,
            zeroline = False
        ),
        yaxis=dict(
            title=dict(text="Observed RTA", font=dict(color="black")),
            range=[0, 1],
            tickmode="array",
            tickvals=[0, 0.5, 1],
            ticktext=["0", "0.5", "1"],
            tickfont=dict(color="black"),
            gridcolor="lightgray",
            linecolor="black"
        ),
        plot_bgcolor="lavender",
        paper_bgcolor="whitesmoke"
    )

    fig = go.Figure(data=[box], layout=layout)
    return fig