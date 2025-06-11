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
import plotly.express as px


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
def plot_rta_boxplot(results_df,
                     box_fill="gray",
                     box_line_color="black",
                     box_width=0.2,
                     jitter=0.5):
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

    # Set up y ticks in a nice way
    max_val = y_vals.max() if len(y_vals) else 1
    ticks2 = max_val*2
    ticks2_int = int(ticks2) + (ticks2 > int(ticks2))
    top_tick = ticks2_int/2

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
            zeroline = False,
            autorange=False,
            fixedrange=True,
            domain=[0.1, 0.9]
        ),
        yaxis=dict(
            title=dict(text="Observed RTA", font=dict(color="black")),
            range=[0, top_tick],
            dtick = 0.5,
            tick0= 0,
            tickfont=dict(color="black"),
            gridcolor="lightgray",
            linecolor="black"
        ),
        plot_bgcolor="#daecf1",
        paper_bgcolor="whitesmoke"
    )

    fig = go.Figure(data=[box], layout=layout)
    return fig

def integrate_metadata(results, metadata):
    """
    Merges the results of RiboMARTA with the available metadata, on same "Run" id.
    :param results: the results from RiboMARTA's run_analysis()
    :param metadata: the metadata, as formatted for use within the Riboseq Data Portal.
    """
    # Creates a slim version of the metadata
    slim_meta = metadata[["Run", "TISSUE", "CELL_LINE", "Strain", "Genotype", "Growth_Condition", "Info", "SampleName", "sample_source", "BioSample"]].copy()
    # assigns the column "Run" as the index
    slim_meta.set_index("Run", inplace=True)

    # merge the two dataframes on index, keep only indexes in results
    merged_df =  pd.merge(slim_meta, results, left_index=True, right_index=True, how="right")

    return merged_df

# Diamond plot function
def diamond_plotter(df, meta="TISSUE"):
    """
    Plots the calculated RTAs by selected meta-information (such as TISSUE or CELL_LINE).
    Significant values are diamonds, while non-significant are circles.
    """
    tissues = df[meta].unique()
    pastel_colors = px.colors.qualitative.Pastel
    tissue_color_map = {t: pastel_colors[i % len(pastel_colors)] for i, t in enumerate(tissues)}

    fig = go.Figure()

    # Circles (significant == False)
    mask_circle = df['significant'] == False
    for tissue in tissues:
        df_t = df[mask_circle & (df[meta] == tissue)]
        fig.add_trace(go.Scatter(
            x=df_t[meta],
            y=df_t['observed_rta'],
            mode='markers',
            marker=dict(
                symbol='circle',
                color=tissue_color_map[tissue],
                size=8,
                line=dict(width=1, color='white')
            ),
            customdata=df_t.index,
            hovertemplate='Tissue: %{x}<br>RTA: %{y}<br>Run: %{customdata}<extra></extra>',
            showlegend=False  # <--- don't show legend here
        ))

    # Diamonds (significant == True)
    mask_diamond = df['significant'] == True
    for tissue in tissues:
        df_t = df[mask_diamond & (df[meta] == tissue)]
        fig.add_trace(go.Scatter(
            x=df_t[meta],
            y=df_t['observed_rta'],
            mode='markers',
            marker=dict(
                symbol='diamond',
                color=tissue_color_map[tissue],
                size=13,
                line=dict(width=1, color='white')
            ),
            customdata=df_t.index,
            hovertemplate='Tissue: %{x}<br>RTA: %{y}<br>Run: %{customdata}<extra></extra>',
            showlegend=False  # <--- don't show legend here
        ))

    # --- Add "fake" traces for legend ---
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='diamond', color='gray', size=13, line=dict(width=1, color='gray')),
        legendgroup='sig',
        showlegend=True,
        name='s. (diamond)'
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='circle', color='gray', size=8, line=dict(width=1, color='gray')),
        legendgroup='nonsig',
        showlegend=True,
        name='n.s. (circle)'
    ))

    # Background, grid, and reference lines
    fig.update_layout(
        margin = dict(b=70),
        plot_bgcolor='#daecf1',
        paper_bgcolor='white',
        title=dict(text = f'Distribution of RTA by {meta}', font = dict(color="black")),
        legend_title_text='Significance:',
        legend = dict(font = dict(color = "black"), title_font = dict(color="black"))
    )

    # axis-level updates:
    fig.update_xaxes(
        automargin=False,
        title_text=meta,
        title_font=dict(color='black'),
        tickfont=dict(color='black'),
        showgrid=True, gridcolor='white',
        tickangle=45
    )
    fig.update_yaxes(
        title_text='Observed RTA',
        title_font=dict(color='black'),
        tickfont=dict(color='black'),
        showgrid=True, gridcolor='white'
    )

    fig.add_hline(y=0, line_color='white', line_dash='solid', line_width=1)

    return fig

# plot to have single transcript plot of the results
def multi_line_plot(df, utr5_length = 0, cds_length = 0, run_option = None, coordinates = None, x=None):
    """
    Plots each column of a dataframe as a line (except for x, if provided).
    Args:
        df: DataFrame (columns to plot as separate lines)
        utr5_length: length of 5' leader
        cds_length: length of cds
        run_option: a list of runs to pass to be plotted
        coordinates: the ranges adopted by the user for the analysis
        x: Optional - Series or column name to use as x-axis (otherwise use index)
        title: Plot title
    """
    palette = px.colors.qualitative.Plotly

    fig = go.Figure()

    if run_option is not None:
        df = df.loc[:, df.columns.isin(run_option)]

    # Decide x values
    if x is not None:
        if isinstance(x, str):
            x_vals = df[x]
        else:
            x_vals = x
        plot_cols = [col for col in df.columns if col != x]
    else:
        x_vals = df.index
        plot_cols = df.columns

    # Add background colours to selected regions
    for i, (region, (start, end)) in enumerate(coordinates.items()):
        fig.add_shape(
            type="rect",
            x0=start,
            x1=end,
            y0=0,
            y1=1,
            yref='paper',  # This makes the rectangle span the full y-axis
            fillcolor= "#FFDFA0",
            line=dict(width=0),
            layer="below"
        )

    # Add labels
    fig.add_annotation(
        x=(coordinates["x_slice"][0] + coordinates["x_slice"][1]) / 2,
        y= 1,
        yref="paper",
        text="1st Region",
        showarrow=False,
        font=dict(size=11, color="black"),
        bgcolor="#FFDFA0",
        borderpad=3
    )
    fig.add_annotation(
        x=(coordinates["y_slice"][0] + coordinates["y_slice"][1]) / 2,
        y= 1,
        yref="paper",
        text="2nd Region",
        showarrow=False,
        font=dict(size=11, color="black"),
        bgcolor="#FFDFA0",
        borderpad=3
    )
    fig.add_annotation(
        x=(coordinates["n_slice"][0] + coordinates["n_slice"][1]) / 2,
        y= 1,
        yref="paper",
        text="Noise Region",
        showarrow=False,
        font=dict(size=11, color="black"),
        bgcolor="#FFDFA0",
        borderpad=3
    )

    # Add a line for each column
    for i, col in enumerate(plot_cols):
        fig.add_trace(go.Bar(
            x=x_vals,
            y=df[col],
            name=col,
            marker=dict(color=palette[i % len(palette)])
            # Optionally: hovertemplate=f"{col}: %{{y}}<extra></extra>"
        ))

    # Add lines to mark CDS
    fig.add_vline(
    x=utr5_length,  # your x-value
    line_dash="dash",  # dashed line
    line_color="black",  # or any color you want
    line_width=2
    )
    fig.add_annotation(
    x=utr5_length,
    y=1,
    yref="paper",
    text="START",
    showarrow=False,
    font=dict(size=12, color="black"),
    bgcolor="#daecf1",
    #bordercolor="black",
    borderpad=4
    )
    # Add lines to mark CDS
    fig.add_vline(
    x=utr5_length+cds_length,  # your x-value
    line_dash="dash",
    line_color="black",
    line_width=2
    )
    fig.add_annotation(
    x=utr5_length+cds_length,
    y=1,
    yref="paper",
    text="STOP",
    showarrow=False,
    font=dict(size=12, color="black"),
    bgcolor="#daecf1",
    #bordercolor="black",
    borderpad=4
    )

    fig.update_layout(
        title="Raw Counts Plot",
        xaxis_title= "Position on Transcript",
        yaxis_title= "Raw counts",
        plot_bgcolor="#daecf1"
    )
    fig.update_yaxes(
        showgrid=True,
        gridcolor="lightgray",
        # ...other y-axis settings...
    )

    #fig.show()
    return fig