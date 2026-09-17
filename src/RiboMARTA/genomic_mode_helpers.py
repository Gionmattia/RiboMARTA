"""
genomic_mode_helpers.py

This module contains functions used by the genomic coordinates mode that were not shared with other modalities.

Author: Gionmattia Carancini
Date: 26-06-2025
"""
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px


# plot to have the three genomic transcript_ranges plotted for the original raw data
def multi_line_plot_genomic(df, run_option = None, coordinates = None, x=None):
    """
    Plots each column of a dataframe as a line (except for x, if provided).
    Args:
        df: DataFrame (columns to plot as separate lines)
        utr5_length: length of 5' leader
        cds_length: length of cds
        run_option: a list of runs to pass to be plotted
        coordinates: the transcript_ranges adopted by the user for the analysis
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
        borderpad=3
    )
    fig.add_annotation(
        x=(coordinates["y_slice"][0] + coordinates["y_slice"][1]) / 2,
        y= 1,
        yref="paper",
        text="2nd Region",
        showarrow=False,
        font=dict(size=11, color="black"),
        borderpad=3
    )
    fig.add_annotation(
        x=(coordinates["n_slice"][0] + coordinates["n_slice"][1]) / 2,
        y= 1,
        yref="paper",
        text="Noise Region",
        showarrow=False,
        font=dict(size=11, color="black"),
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

    for key in coordinates.keys():
        # Add lines to mark CDS
        fig.add_vline(
            x=coordinates[key][0],  # your x-value
            line_dash="solid",  # dashed line
            line_color="black",  # or any color you want
            line_width=2
        )

        fig.add_vline(
            x=coordinates[key][1],  # your x-value
            line_dash="solid",  # dashed line
            line_color="black",  # or any color you want
            line_width=2
        )


    fig.update_layout(
        title="Raw Counts Plot",
        xaxis_title= "Position over the three genomic transcript_ranges",
        yaxis_title= "Raw counts",
        plot_bgcolor="#daecf1"
    )
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(
        showgrid=True,
        gridcolor="lightgray"
    )

    #fig.show()
    return fig


# Functions to make the Aggregated Plot by Metadata

# Helper
def _extract_runs_by_meta(results_meta, agg_meta_option):
    """
    Internal helper, used to retrieve all runs divided by selected agg_meta_option

    :param results_meta: Datafarme with runs as indexes and metadata
    :param agg_meta_option: user-specified paramter rto group (and later plot) runs by
    :return: dictionary where each key is a value of the agg_meta_option selected, value is the list of runs.
    """
    grouped_runs = {}

    for group_name, group_df in results_meta.groupby(agg_meta_option):
        grouped_runs[group_name] = group_df.index.tolist()

    return grouped_runs


# Function to generate aggregated (by metadata) profiles.
def aggregator(input_df, results_meta, agg_meta_option):
    """
    Normalises the raw data and aggregates together studies based on metadata field given as input

    :param input_df: The raw data
    :param results_meta: Datafarme with runs results and metadata
    :param agg_meta_option: the metadata field to use in the aggregation
    :return:
    """

    grouped_runs = _extract_runs_by_meta(results_meta, agg_meta_option)
    # normalize by total counts per run
    df = input_df / input_df.sum()

    # Set up output
    temp_dict = {}

    # Aggregate runs based on the grouped_runs dictionary
    for meta, run_list in grouped_runs.items():
        filtered_df_cols = df.loc[:, run_list]
        group_sum_series = filtered_df_cols.sum(axis=1)  # summing row-wise

        temp_dict[meta] = group_sum_series

    # Return dataframe of global aggregates based on specified meta
    agg_output_df = pd.DataFrame(temp_dict)
    agg_output_df["Position"] = input_df["Position"]

    return agg_output_df


# plot aggregated coverage for selected meta conditions
def genom_agg_line_plot(agg_output_df, aggregations = None, coordinates = None, x=None):
    """
    Plots each column of a dataframe as a line (except for x, if provided).
    Args:
        agg_output_df: DataFrame result of aggregator
        utr5_length: length of 5' leader
        cds_length: length of cds
        aggregations: a list of metadata options for the field runs have been grouped by aggregator()
        coordinates: the transcript_ranges adopted by the user for the analysis
        x: Optional - Series or column name to use as x-axis (otherwise use index)
    """

    # First, normalise the data for comparability

    # Second, extract (as dictionary) the studies matching aggregations. the value is the key, the value is the list of studies.

    # For key in dict, aggregate the studies matching the list from df_normalised and append the result as a new column to "df"

    df = agg_output_df.loc[:,aggregations].copy()

    # Create a new dataframe, where we have a column for each different aggregation


    palette = px.colors.qualitative.Plotly

    fig = go.Figure()

    if aggregations is not None:
        df = df.loc[:, df.columns.isin(aggregations)]

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

    df["Position"] = agg_output_df["Position"]

    # Add labels
    fig.add_annotation(
        x=(coordinates["x_slice"][0] + coordinates["x_slice"][1]) / 2,
        y= 1,
        yref="paper",
        text="1st Region",
        showarrow=False,
        font=dict(size=11, color="black"),
        borderpad=3
    )
    fig.add_annotation(
        x=(coordinates["y_slice"][0] + coordinates["y_slice"][1]) / 2,
        y= 1,
        yref="paper",
        text="2nd Region",
        showarrow=False,
        font=dict(size=11, color="black"),
        borderpad=3
    )
    fig.add_annotation(
        x=(coordinates["n_slice"][0] + coordinates["n_slice"][1]) / 2,
        y= 1,
        yref="paper",
        text="Noise Region",
        showarrow=False,
        font=dict(size=11, color="black"),
        borderpad=3
    )

    # Add a line for each column
    for i, col in enumerate(plot_cols):
        fig.add_trace(go.Bar(
            x=x_vals,
            y=df[col],
            name=col,
            marker=dict(color=palette[i % len(palette)]),
            customdata = df["Position"],
            hovertemplate=("Position: %{customdata}<br>")
        ))

    for key in coordinates.keys():
        # Add lines to mark CDS
        fig.add_vline(
            x=coordinates[key][0],  # your x-value
            line_dash="solid",  # dashed line
            line_color="black",  # or any color you want
            line_width=2
        )

        fig.add_vline(
            x=coordinates[key][1],  # your x-value
            line_dash="solid",  # dashed line
            line_color="black",  # or any color you want
            line_width=2
        )

    fig.update_layout(
        title="Aggregated Counts Plot",
        xaxis_title= "Position over the three genomic transcript_ranges",
        yaxis_title= "Raw counts",
        plot_bgcolor="#daecf1"
    )
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(
        showgrid=True,
        gridcolor="lightgray",
    )

    #fig.show()
    return fig