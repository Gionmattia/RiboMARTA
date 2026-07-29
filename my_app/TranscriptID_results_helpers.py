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
        min_y_range, max_y_range = 0, max(1, results_df[y_column].max() + 0.05)  # Should use a np.max and select either 1.15 or the highest among the data)

    filtered = results_df[y_column].dropna()
    y_vals = filtered.values

    # Safety check (against cases where there are no results)
    if len(y_vals) == 0:
        fig = go.Figure()
        fig.update_layout(
            title=dict(text=plot_title, font=dict(color="black")),
            plot_bgcolor="white",
            paper_bgcolor="white",
            annotations=[dict(
                x=0.5, y=0.5,
                xref="paper", yref="paper",
                text="No data available for this selection",
                showarrow=False,
                font=dict(color="gray", size=14),
            )]
        )
        return fig

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


# Function to make diamond plot with CI whiskers by selected metadata field
def new_diamond_plotter(df, meta="TISSUE", logs=False, jitter_amount=0.4, show_nonsignificant=False):
    """
    Plots the calculated RTAs by selected meta-information (such as TISSUE or CELL_LINE).
    Significant values are diamonds, while non-significant are circles (optional).
    Each point has confidence interval whiskers.
    Points are jittered horizontally within each tissue group.

    Args:
        df: DataFrame with RTA results
        meta: Metadata column to group by (default: "TISSUE")
        logs: If True, plot log2RTA values (default: False)
        jitter_amount: Amount of horizontal jitter within groups (default: 0.4)
        show_nonsignificant: If True, show non-significant points as circles (default: False)
    """
    df = df.copy()
    if not show_nonsignificant:
        df = df[df['significant'] == True]

    tissues = df[meta].unique()
    pastel_colors = px.colors.qualitative.Pastel
    tissue_color_map = {t: pastel_colors[i % len(pastel_colors)] for i, t in enumerate(tissues)}

    if logs == True:
        y_column = "log2RTA"
        df["log2RTA_ci_lower_95%"] = np.log2(df["RTA_ci_lower_95%"])
        df["log2RTA_ci_upper_95%"] = np.log2(df["RTA_ci_upper_95%"])
        ci_lower = "log2RTA_ci_lower_95%"
        ci_upper = "log2RTA_ci_upper_95%"
        plot_title = f"Diamond plot of log2 RTAs across datasets, by {meta.lower()}"
        tot_counts = "Tot Counts"
    else:
        y_column = "RTA"
        ci_lower = "RTA_ci_lower_95%"
        ci_upper = "RTA_ci_upper_95%"
        plot_title = f"Diamond plot of RTAs across datasets, by {meta.lower()}"
        tot_counts = "Tot Counts"

    fig = go.Figure()

    # Create a mapping of tissue to numeric position
    tissue_positions = {tissue: i for i, tissue in enumerate(tissues)}

    # Add jittered positions to dataframe
    df = df.copy()
    df['x_pos'] = df[meta].map(tissue_positions)

    # Add random jitter within each tissue group
    np.random.seed(42)  # For reproducibility
    for tissue in tissues:
        mask = df[meta] == tissue
        n_points = mask.sum()
        jitter = np.random.uniform(-jitter_amount, jitter_amount, n_points)
        df.loc[mask, 'x_jitter'] = df.loc[mask, 'x_pos'] + jitter

    # Add error bars first (so they appear behind the points)
    for tissue in tissues:
        df_tissue = df[df[meta] == tissue]

        # Add error bars for this tissue
        for idx, row in df_tissue.iterrows():
            fig.add_trace(go.Scatter(
                x=[row['x_jitter'], row['x_jitter']],
                y=[row[ci_lower], row[ci_upper]],
                mode='lines',
                line=dict(color='gray', width=1.5),
                showlegend=False,
                hoverinfo='skip'
            ))

    # Circles (significant == False). By default hidden.
    mask_circle = df['significant'] == False
    for tissue in tissues:
        df_t = df[mask_circle & (df[meta] == tissue)]
        if len(df_t) > 0:
            fig.add_trace(go.Scatter(
                x=df_t['x_jitter'],
                y=df_t[y_column],
                mode='markers',
                marker=dict(
                    symbol='circle',
                    color=tissue_color_map[tissue],
                    size=8,
                    line=dict(width=1, color='white')
                ),
                customdata=np.column_stack((df_t[meta], # customdata[0]
                                            df_t.index, # customdata[1]
                                            df_t[ci_lower].round(2), # customdata[2]
                                            df_t[ci_upper].round(2), # customdata[3]
                                            df_t[tot_counts])),  # customdata[4]

                hovertemplate='Tissue: %{customdata[0]}<br>'
                              'RTA: %{y:.2f}<br>'
                              'CI: %{customdata[2]} - %{customdata[3]}<br>'
                              'Run: %{customdata[1]}<br>'
                              'Tot Counts: %{customdata[4]}<extra></extra>',
                showlegend=False
            ))

    # Diamonds (significant == True)
    mask_diamond = df['significant'] == True
    for tissue in tissues:
        df_t = df[mask_diamond & (df[meta] == tissue)]
        fig.add_trace(go.Scatter(
            x=df_t['x_jitter'],
            y=df_t[y_column],
            mode='markers',
            marker=dict(
                symbol='diamond',
                color=tissue_color_map[tissue],
                size=13,
                line=dict(width=1, color='white')
            ),
            customdata=np.column_stack((df_t[meta],  # customdata[0]
                                        df_t.index,  # customdata[1]
                                        df_t[ci_lower].round(2),  # customdata[2]
                                        df_t[ci_upper].round(2),  # customdata[3]
                                        df_t[tot_counts])),  # customdata[4]
            hovertemplate='Tissue: %{customdata[0]}<br>'
                          'RTA: %{y:.2f}<br>'
                          'CI: %{customdata[2]} - %{customdata[3]}<br>'
                          'Run: %{customdata[1]}<br>'
                          'Tot Counts: %{customdata[4]}<extra></extra>',
            showlegend=False
        ))

    # --- Add "fake" traces for legend ---
    if show_nonsignificant:
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
        margin=dict(b=70),
        plot_bgcolor='#daecf1',
        paper_bgcolor='white',
        title=dict(text=plot_title, font=dict(color="black")),
        legend_title_text='Significance:',
        legend=dict(font=dict(color="black"), title_font=dict(color="black")),
        width=max(600, len(tissues) * 100)
    )

    # axis-level updates with categorical labels
    fig.update_xaxes(
        automargin=False,
        title_text=meta.lower(),
        title_font=dict(color='black'),
        tickfont=dict(color='black'),
        showgrid=True, gridcolor='white',
        tickangle=45,
        tickmode='array',
        tickvals=list(tissue_positions.values()),
        ticktext=list(tissue_positions.keys())
    )
    fig.update_yaxes(
        title_text=y_column,
        title_font=dict(color='black'),
        tickfont=dict(color='black'),
        showgrid=True, gridcolor='white'
    )

    fig.add_hline(y=0, line_color='white', line_dash='solid', line_width=1)

    return fig


import numpy as np
import pandas as pd
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool, FixedTicker, Span
from bokeh.models import Legend

def new_diamond_plotter_bokeh(df, meta="TISSUE", logs=False, jitter_amount=0.4,
                              show_nonsignificant=False, show_n_members = False):
    """
    Plots the calculated RTAs by selected meta-information (such as TISSUE or CELL_LINE).
    Significant values are diamonds, while non-significant are circles.
    Each point estimate has its CI wiskers.
    Points are jittered horizontally within each group.

    Args:
        df: DataFrame with RTA results.
        meta: Metadata column to group by (default: "TISSUE").
        logs: If True, plot log2RTA values (default: False). Expects a "log2RTA"
              column to already exist in df (matching the original).
        jitter_amount: Amount of horizontal jitter within groups (default: 0.4).
        show_nonsignificant: If True, show non-significant points as circles
                             (default: False).

    Returns:
        A Bokeh figure object.
    """
    df = df.copy()
    if not show_nonsignificant:
        df = df[df['significant'] == True]

    tissues = df[meta].unique()

    # Color Palette (added manually)
    pastel_colors = [
        'rgb(102, 197, 204)', 'rgb(246, 207, 113)', 'rgb(248, 156, 116)',
        'rgb(220, 176, 242)', 'rgb(135, 197, 95)', 'rgb(158, 185, 243)',
        'rgb(254, 136, 177)', 'rgb(201, 219, 116)', 'rgb(139, 224, 164)',
        'rgb(180, 151, 231)', 'rgb(179, 179, 179)'
    ]
    tissue_color_map = {t: pastel_colors[i % len(pastel_colors)]
                        for i, t in enumerate(tissues)}


# This might not work if CI are different than 95%... well it stays for now
    if logs:
        y_column = "log2RTA"
        df["log2RTA_ci_lower_95%"] = np.log2(df["RTA_ci_lower_95%"])
        df["log2RTA_ci_upper_95%"] = np.log2(df["RTA_ci_upper_95%"])
        ci_lower = "log2RTA_ci_lower_95%"
        ci_upper = "log2RTA_ci_upper_95%"
        plot_title = f"Diamond plot of log2 RTAs across datasets, by {meta.lower()}"
    else:
        y_column = "RTA"
        ci_lower = "RTA_ci_lower_95%"
        ci_upper = "RTA_ci_upper_95%"
        plot_title = f"Diamond plot of RTAs across datasets, by {meta.lower()}"

    tot_counts = "Tot Counts"

    # Map each tissue to a numeric x position (categorical axis emulated numerically).
    tissue_positions = {tissue: i for i, tissue in enumerate(tissues)}
    df['x_pos'] = df[meta].map(tissue_positions)

    # Add random jitter within each tissue group.
    # Same seed + same iteration order as the original -> identical jitter values.
    np.random.seed(42)
    df['x_jitter'] = np.nan
    for tissue in tissues:
        mask = df[meta] == tissue
        n_points = int(mask.sum())
        jitter = np.random.uniform(-jitter_amount, jitter_amount, n_points)
        df.loc[mask, 'x_jitter'] = df.loc[mask, 'x_pos'] + jitter

    df['color'] = df[meta].map(tissue_color_map)

    p = figure(
        frame_height=450,
        sizing_mode="stretch_width",
        title=plot_title,
        x_axis_label=meta.lower(),
        y_axis_label=y_column,
        tools="pan,box_zoom,wheel_zoom,reset,save",
        background_fill_color='#daecf1',
    )

    # --- CI whiskers first (so they sit behind the markers) ---
    whisker_source = ColumnDataSource(dict(
        x=df['x_jitter'].to_numpy(),
        y_lower=df[ci_lower].to_numpy(),
        y_upper=df[ci_upper].to_numpy(),
    ))
    p.segment('x', 'y_lower', 'x', 'y_upper', source=whisker_source,
              line_color='gray', line_width=1.5)
    if show_n_members == False:
        def make_source(sub):
            return ColumnDataSource(dict(
                x=sub['x_jitter'].to_numpy(),
                y=sub[y_column].to_numpy(),
                color=sub['color'].to_numpy(),
                meta=sub[meta].to_numpy(),
                run=sub.index.to_numpy(),
                ci_lower=sub[ci_lower].round(2).to_numpy(),
                ci_upper=sub[ci_upper].round(2).to_numpy(),
                tot_counts=sub[tot_counts].to_numpy()
            ))
    elif show_n_members == True:
        def make_source(sub):
            return ColumnDataSource(dict(
                x=sub['x_jitter'].to_numpy(),
                y=sub[y_column].to_numpy(),
                color=sub['color'].to_numpy(),
                meta=sub[meta].to_numpy(),
                run=sub.index.to_numpy(),
                ci_lower=sub[ci_lower].round(2).to_numpy(),
                ci_upper=sub[ci_upper].round(2).to_numpy(),
                tot_counts=sub[tot_counts].to_numpy(),
                count = sub["count"].to_numpy()
            ))

    hover_renderers = []

    # --- Circles (significant == False); only present when show_nonsignificant ---
    df_circle = df[df['significant'] == False]
    if len(df_circle) > 0:
        circle_source = make_source(df_circle)
        circle_r = p.scatter('x', 'y', source=circle_source, marker='circle',
                             size=8, fill_color='color', line_color='white',
                             line_width=1)
        hover_renderers.append(circle_r)

    # --- Diamonds (significant == True) ---
    df_diamond = df[df['significant'] == True]
    diamond_source = make_source(df_diamond)
    diamond_r = p.scatter('x', 'y', source=diamond_source, marker='diamond',
                          size=25, fill_color='color', line_color='white',
                          line_width=1)
    hover_renderers.append(diamond_r)

    # --- Hover (mirrors the original hovertemplate fields) ---

    if show_n_members == False:
        hover = HoverTool(renderers=hover_renderers, tooltips=[
            (meta, "@meta"),
            ("RTA", "@y{0.00}"),
            ("CI", "@ci_lower - @ci_upper"),
            ("Run", "@run"),
            ("Tot Counts", "@tot_counts"),
        ])
    elif show_n_members == True:
        hover = HoverTool(renderers=hover_renderers, tooltips=[
            (meta, "@meta"),
            ("RTA", "@y{0.00}"),
            ("CI", "@ci_lower - @ci_upper"),
            ("Run", "@run"),
            ("Tot Counts", "@tot_counts"),
            ("N Datasets", "@count"),
        ])

    p.add_tools(hover)

    # --- Legend (significance), only when non-significant points are shown ---
    if show_nonsignificant:
        # Empty renderers: they contribute a legend entry without drawing points,
        # matching the original's x=[None], y=[None] "fake" traces.
        d = p.scatter(x=[float('nan')], y=[float('nan')], marker='diamond', size=25, fill_color='black',
                      line_color='black')
        c = p.scatter(x=[float('nan')], y=[float('nan')], marker='circle', size=8, fill_color='black',
                      line_color='black')
        legend = Legend(items=[('s.', [d]), ('n.s.', [c])], title='Significance:')
        legend.label_text_color = 'black'
        legend.title_text_color = 'black'

        p.add_layout(legend, "right")

    # --- y = 0 reference line (white, solid) ---
    p.add_layout(Span(location=0, dimension='width',
                      line_color='white', line_width=1))

    # --- Categorical tick labels on the numeric x-axis ---
    p.xaxis.ticker = FixedTicker(ticks=list(tissue_positions.values()))
    p.xaxis.major_label_overrides = {i: str(t) for t, i in tissue_positions.items()}
    p.xaxis.major_label_orientation = np.pi / 4  # 45 degrees

    # --- Cosmetics to match the original (white grid on pale-blue panel) ---
    p.xgrid.grid_line_color = 'white'
    p.ygrid.grid_line_color = 'white'
    p.title.text_color = 'black'
    p.xaxis.axis_label_text_color = 'black'
    p.yaxis.axis_label_text_color = 'black'
    p.xaxis.major_label_text_color = 'black'
    p.yaxis.major_label_text_color = 'black'

    return p

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
        text="Baseline",
        showarrow=False,
        font=dict(size=11, color="black"),
        bgcolor="#FFDFA0",
        borderpad=3
    )
    fig.add_annotation(
        x=(coordinates["y_slice"][0] + coordinates["y_slice"][1]) / 2,
        y= 1,
        yref="paper",
        text="Target",
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
        title=f"Raw Counts Plot for {str(run_option).replace('[','').replace(']','')}",
        xaxis_title= "Position on Transcript",
        yaxis_title= "Raw counts",
        plot_bgcolor="#daecf1"
    )
    fig.update_yaxes(
        showgrid=True,
        gridcolor="lightgray",
    )

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
        title=f"Aggregated Counts for {str(aggregations).replace('[','').replace(']','')}",
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