import pandas as pd
from scipy.stats import chi2_contingency
import numpy as np
import plotly.graph_objects as go
from statsmodels.stats.multitest import multipletests

# Does the breakdown (across UTRs and CDS) from the raw counts. Selects only runs that made through the RiboMARTA analysis
def get_breakdowns(all_studies_df, txid_results, utr_5, cds):
    good_runs = txid_results.index.to_list()
    df = all_studies_df[good_runs].copy()
    # Compute the three totals by region...
    leader_totals = df.iloc[:utr_5].sum(axis=0)
    leader_totals.name = "UTR_5"

    cds_totals = df.iloc[utr_5:utr_5 + cds].sum(axis=0)
    cds_totals.name = "CDS"

    trailer_totals = df.iloc[utr_5 + cds:].sum(axis=0)
    trailer_totals.name = "UTR_3"

    breakdown = pd.concat([leader_totals, cds_totals, trailer_totals], axis=1)

    breakdown["Run"] = breakdown.index

    df = breakdown.drop("Run", axis=1)

    return df

# Does Chi-squared
def execute_chi_squared(contingency_table):
    # Chi squared
    chi2, p_value, dof, expected_frequencies = chi2_contingency(contingency_table)
    # Cramer's V correction
    n = contingency_table.sum().sum()  # Total sample size
    r, c = contingency_table.shape  # Number of rows and columns

    # Calculate the minimum of (r-1) and (c-1)
    min_dim_minus_1 = min(r - 1, c - 1)

    # Handle cases where min_dim_minus_1 could be zero (e.g., 1xN or Nx1 table)
    if min_dim_minus_1 == 0:
        cramers_v = 0.0  # If one dimension is 1, there's no variability for association
    else:
        cramers_v = np.sqrt(chi2 / (n * min_dim_minus_1))

    std_residuals = (contingency_table - expected_frequencies) / np.sqrt(expected_frequencies)

    return chi2, p_value, dof, cramers_v, std_residuals

# Allows visualisation of std_residuals. USeful to see what differs the most.
def plot_chisquared_residuals(std_residuals):
    fig = go.Figure(
        data=go.Heatmap(
            z=std_residuals.values,
            x=std_residuals.columns,
            y=std_residuals.index,
            colorscale="RdBu",
            zmin=-abs(std_residuals.values).max(),
            zmax=abs(std_residuals.values).max(),
            colorbar=dict(title="Residual"),
            text=std_residuals.round(2).astype(str),
            hovertemplate="Row: %{y}<br>Col: %{x}<br>Residual: %{z:.2f}<extra></extra>"
        )
    )
    fig.update_layout(
        title="Standardized Residuals Heatmap",
        xaxis_title="Category",
        yaxis_title="Group"
    )

    return fig


def get_breakdowns_slices(all_studies_df, txid_results, coordinates):
    good_runs = txid_results.index.to_list()
    df = all_studies_df[good_runs].copy()

    # Compute the three totals by region...

    first_region = df.iloc[coordinates["x_slice"][0]:coordinates["x_slice"][1]].sum(axis=0)
    first_region.name = "First Region"

    second_region = df.iloc[coordinates["y_slice"][0]:coordinates["y_slice"][1]].sum(axis=0)
    second_region.name = "Second Region"

    breakdown = pd.concat([first_region, second_region], axis=1)

    breakdown["Run"] = breakdown.index

    df = breakdown.drop("Run", axis=1)

    return df

def get_breakdowns_percent(all_studies_df, txid_results, utr_5, cds):
    good_runs = txid_results.index.to_list()
    df = all_studies_df[good_runs].copy()
    # Divide each row-column by column total
    df = df / df.sum(axis=0)

    # Compute the three totals by region...
    leader_totals = df.iloc[:utr_5].sum(axis=0)
    leader_totals = leader_totals * 100
    leader_totals.name = "UTR_5"

    cds_totals = df.iloc[utr_5:utr_5 + cds].sum(axis=0)
    cds_totals = cds_totals * 100
    cds_totals.name = "CDS"

    trailer_totals = df.iloc[utr_5 + cds:].sum(axis=0)
    trailer_totals = trailer_totals * 100
    trailer_totals.name = "UTR_3"

    breakdown = pd.concat([leader_totals, cds_totals, trailer_totals], axis=1)

    breakdown["Run"] = breakdown.index

    df = breakdown.drop("Run", axis=1)

    return df


def get_breakdowns_percent_with_GA(all_studies_df, txid_results, utr_5, cds):
    good_runs = txid_results.index.to_list()
    df = all_studies_df[good_runs].copy()
    # Divide each row-column by column total
    df = df / df.sum(axis=0)
    # Compute normalised global aggregate
    df["Global Aggregate"] = df.sum(axis = 1)
    df["Global Aggregate"] = df["Global Aggregate"]/df["Global Aggregate"].sum(axis = 0)
    # Compute the three totals by region...
    leader_totals = df.iloc[:utr_5].sum(axis=0)
    leader_totals = leader_totals * 100
    leader_totals.name = "UTR_5"

    cds_totals = df.iloc[utr_5:utr_5 + cds].sum(axis=0)
    cds_totals = cds_totals * 100
    cds_totals.name = "CDS"

    trailer_totals = df.iloc[utr_5 + cds:].sum(axis=0)
    trailer_totals = trailer_totals * 100
    trailer_totals.name = "UTR_3"

    breakdown = pd.concat([leader_totals, cds_totals, trailer_totals], axis=1)

    breakdown["Run"] = breakdown.index

    df = breakdown.drop("Run", axis=1)

    return df


def perform_pairwise_chi(df2, baseline_idx):
    """
    Perform pairwise Chi-squared tests of all rows in df2 versus the baseline row.
    Returns a DataFrame with chi2, p-value, dof, and corrected p-values.
    """
    # Ensure only numeric columns
    df2_num = df2.select_dtypes(include='number')
    results = []

    for idx, row in df2_num.iterrows():
        if idx == baseline_idx:
            continue
        table = pd.DataFrame([df2_num.loc[baseline_idx], row])
        chi2, p, dof, expected = chi2_contingency(table)
        results.append({
            "row": idx,
            "chi2": chi2,
            "p-value": p,
            "dof": dof
        })

    results_df = pd.DataFrame(results)

    if not results_df.empty:
        # Bonferroni correction
        bonf_corrected = multipletests(results_df['p-value'], method='bonferroni')[1]
        # Benjamini-Hochberg FDR correction
        fdr_corrected = multipletests(results_df['p-value'], method='fdr_bh')[1]
        results_df['bonferroni_p'] = bonf_corrected
        results_df['fdr_p'] = fdr_corrected

    return results_df


def pairwise_comparison_plotter(results_df):
    # Make sure 'row' is string for the axis
    x_labels = results_df['row'].astype(str).values
    y_labels = ['bonferroni_p']  # Only one row, since it's a 1D heatmap
    z_values = [results_df['bonferroni_p'].values]  # 2D array: 1 x N

    # Create custom hover text
    hover_text = [
        [f"Row: {x}<br>bonferroni_p: {p:.4g}" for x, p in zip(x_labels, z_values[0])]
    ]

    fig = go.Figure(data=go.Heatmap(
        z=z_values,
        x=x_labels,
        y=y_labels,
        colorscale='Viridis',
        colorbar=dict(title='bonferroni_p'),
        zmin=0,
        zmax=1,
        hoverinfo='text',
        text=hover_text
    ))

    fig.update_layout(
        title="Pairwise Chi-squared bonferroni_p-values (vs baseline)",
        xaxis_title="Compared Row",
        yaxis_visible=False,  # Hide y axis label for a single row
        yaxis_showticklabels=False,
        margin=dict(l=40, r=40, t=50, b=50),
        width=60*len(x_labels)+120,
        height=200,
    )

    return fig


def better_plotter(df):
    region_cols = ['UTR_5', 'CDS', 'UTR_3']
    x_labels = region_cols
    y_labels = df.index.astype(str)
    z = df[region_cols].values  # n_rows x 3

    # Hover text: include both region value and p-value for the row
    hovertext = np.array([
        [f"{row}<br>{col}: {val:.2f}% of Translational Activity <br>Adjusted (Bonf) p-value: {pval:.3g}"
         for col, val in zip(region_cols, row_vals)]
        for row, row_vals, pval in zip(y_labels, z, df['bonferroni_p'])
    ])

    # Find significant rows (star marker)
    star_rows = df.index[df['bonferroni_p'] < 0.05].tolist()
    star_y = [y_labels.tolist().index(str(idx)) for idx in star_rows]

    fig = go.Figure()

    # Heatmap
    fig.add_trace(go.Heatmap(
        z=z,
        x=x_labels,
        y=y_labels,
        colorscale='Viridis',
        zmin=0, zmax=100,
        hoverinfo='text',
        text=hovertext,
        showscale=True,
        colorbar=dict(
            title='% Region',
            thickness = 12,
            len = 0.1,
            y = 0.9)
    ))

    # Stars for significant p-values
    if star_y:
        fig.add_trace(go.Scatter(
            x=[x_labels[-1]]*len(star_y),  # Place at rightmost column
            y=[y_labels[i] for i in star_y],
            mode='text',
            text=['⭐']*len(star_y),
            textfont=dict(color='red', size=20),
            showlegend=False,
            hoverinfo='skip'
        ))

    fig.update_layout(
        title="Heatmap of Regional Translational Acitivity (%)<br>",
        xaxis_title="Regional activity",
        yaxis_title="Run",
        width=420,
        height=30 * len(df),
        margin=dict(l=60, r=40, t=60, b=40)
    )

    fig.update_xaxes(side='top')

    return fig