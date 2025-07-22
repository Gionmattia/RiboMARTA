import pandas as pd
from scipy.stats import chi2_contingency
import numpy as np
import plotly.graph_objects as go

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

