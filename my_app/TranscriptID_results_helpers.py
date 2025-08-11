"""
TranscriptID_results_helpers.py

This module contains all the functions used to plot and explore the results from YTranscriptID mode.
Author: Gionmattia Carancini
Date: 16-06-2025
"""
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import numpy as np

# Function to make the boxplot of RTAs
def plot_rta_boxplot(results_df,
                     logs=False,
                     box_fill="#00bf7d",
                     box_line_color="gray",
                     box_width=0.2,
                     jitter=0.5):
    """
    Interactive box+points plot of “RTA” (where “significant” == True),
    with hover labels, custom box color, and controlled box width.
    """
    #
    if logs == True:
        y_column = "log2RTA"
        plot_title = "Boxplot of log2 RTAs across datasets"
        min_y_range, max_y_range = results_df[y_column].min() - 0.25, results_df[y_column].max() + 0.25

    else:
        y_column = "RTA"
        plot_title = "Boxplot of observed RTAs across datasets"
        min_y_range, max_y_range = 0, 1

    filtered = results_df[y_column].dropna()
    y_vals = filtered.values
    labels = filtered.index.astype(str).tolist()
    x_vals = [1] * len(y_vals)

    # Set up y ticks in a nice way
    max_val = y_vals.max() if len(y_vals) else 1
    ticks2 = max_val*2
    ticks2_int = int(ticks2) + (ticks2 > int(ticks2))
    # top_tick = ticks2_int/2

    # Calculate descriptive statistics and text
    median_val = np.median(y_vals) if len(y_vals) else np.nan
    mean_val = np.mean(y_vals) if len(y_vals) else np.nan
    n_val = len(y_vals)
    iqr = np.percentile(y_vals, 75) - np.percentile(y_vals, 25)
    stats_text = (
        f"<b>Descriptive Stats</b><br>"
        f"n = {n_val}<br>"
        f"Mean = {mean_val:.2f}<br>"
        f"Median = {median_val:.2f}<br>"
        f"IQR = {iqr:.2f}"
    )

    box = go.Box(
        x=x_vals,
        y=y_vals,
        boxpoints="all",
        jitter=jitter,
        pointpos=0,
        width=box_width,
        fillcolor= box_fill,
        line=dict(color=box_line_color),
        marker=dict(color="black", opacity=0.5),
        text=labels,
        hoverinfo="text+y",
        name="",
        showlegend=False
    )

    layout = go.Layout(
        title=dict(
            text= plot_title,
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
            title=dict(text=y_column, font=dict(color="black")),
            range=[min_y_range, max_y_range],
            dtick = 0.5,
            tick0= 0,
            tickfont=dict(color="black"),
            gridcolor="lightgray",
            linecolor="black"
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        annotations=[
            # Code for median label
            dict(
                x=0.90, y=0.97,
                xref="paper",
                yref="paper",
                text=stats_text,
                showarrow=False,
                align = "left",
                #ax=40,
                #ay=-30,
                font=dict(color="black", size=12),
                bgcolor="white",
                bordercolor="black",
                borderpad=6,

            )
        ]
    )

    fig = go.Figure(data=[box], layout=layout)
    return fig


# Function to make diamond plot my selected metadata field (defo is TISSUE)
def diamond_plotter(df, meta="TISSUE", logs = False):
    """
    Plots the calculated RTAs by selected meta-information (such as TISSUE or CELL_LINE).
    Significant values are diamonds, while non-significant are circles.
    """
    tissues = df[meta].unique()
    pastel_colors = px.colors.qualitative.Pastel
    tissue_color_map = {t: pastel_colors[i % len(pastel_colors)] for i, t in enumerate(tissues)}

    if logs == True:
        y_column = "log2RTA"
        plot_title = f"Diamond plot of log2 RTAs across datasets, by {meta.lower()}"
        #min_y_range, max_y_range = results_df[y_column].min() - 0.25, results_df[y_column].max() + 0.25

    else:
        y_column = "RTA"
        plot_title = f"Diamond plot of RTAs across datasets, by {meta.lower()}"
        #min_y_range, max_y_range = 0, 1


    fig = go.Figure()

    # Circles (significant == False)
    mask_circle = df['significant'] == False
    for tissue in tissues:
        df_t = df[mask_circle & (df[meta] == tissue)]
        fig.add_trace(go.Scatter(
            x=df_t[meta],
            y=df_t[y_column],
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
            y=df_t[y_column],
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
        name='s.'
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='circle', color='gray', size=8, line=dict(width=1, color='gray')),
        legendgroup='nonsig',
        showlegend=True,
        name='n.s.'
    ))

    # Background, grid, and reference lines
    fig.update_layout(
        margin = dict(b=70),
        plot_bgcolor='#daecf1',
        paper_bgcolor='white',
        title=dict(text = plot_title, font = dict(color="black")),
        legend_title_text='Significance:',
        legend = dict(font = dict(color = "black"), title_font = dict(color="black"))
    )

    # axis-level updates:
    fig.update_xaxes(
        automargin=False,
        title_text=meta.lower(),
        title_font=dict(color='black'),
        tickfont=dict(color='black'),
        showgrid=True, gridcolor='white',
        tickangle=45
    )
    fig.update_yaxes(
        title_text= y_column,
        title_font=dict(color='black'),
        tickfont=dict(color='black'),
        showgrid=True, gridcolor='white'
    )

    fig.add_hline(y=0, line_color='white', line_dash='solid', line_width=1)

    return fig


# plot to have single transcript plot of the original raw data
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

    return agg_output_df


# plot aggregated coverage for selected meta conditions
def agg_line_plot(agg_output_df,  utr5_length = 0, cds_length = 0, aggregations = None, coordinates = None, x=None):
    """
    Plots each column of a dataframe as a line (except for x, if provided).
    Args:
        agg_output_df: DataFrame result of aggregator
        utr5_length: length of 5' leader
        cds_length: length of cds
        aggregations: a list of metadata options for the field runs have been grouped by aggregator()
        coordinates: the ranges adopted by the user for the analysis
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
        title="Aggregated Counts Plot",
        xaxis_title= "Position on Transcript",
        yaxis_title= "Normalised counts (aggregated by tissue)",
        plot_bgcolor="#daecf1"
    )
    fig.update_yaxes(
        showgrid=True,
        gridcolor="lightgray",
        # ...other y-axis settings...
    )

    #fig.show()
    return fig