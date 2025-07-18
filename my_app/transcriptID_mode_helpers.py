# transcriptID_mode_helpers.py

import streamlit as st
import pyBigWig
import pyranges as pr
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing as mp
mp.set_start_method("fork", force=True)
import plotly.graph_objects as go
import random
from module0 import _extract_one
import time
from Bio.Seq import Seq
from bokeh.models import Span, Label
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import HoverTool
from bokeh.models import PanTool

# Function to load the ranges
def extract_tx_pyranges(gr, txid):
    """
    Function to load the gtf annotation and extract the ranges for selected transcript

    :param gr: loaded gtf, as pyranges object
    :param txid: the stable (ensebml) transcript id
    :return: a dataframe with ranges for each feature of txid
    """
    # Read the gtf file
    transcript_ranges = gr[gr.transcript_id == txid]
    if len(transcript_ranges) == 0:
        raise KeyError(f"transcript id {txid} was not found in the annotation file.")

    return transcript_ranges


# Function to format the transcript ranges, for extraction of coverage from BigWigs
def regions_dataframe_formatter(transcript_ranges, return_individual_regions = False):
    """
    :param transcript_ranges: a pyranges.pyranges_main.PyRanges object
    :param return_individual_regions: whether to return individual regions (3', 5', cds) dataframe coordinates or as a single one
    :return: a singe dataframe, with the coordinates for each exon fragment making up CDS and utrs or the dataframes with coordinates for utr5, CDS, and utr3 separately.

    Regarding this aspect: On the + strand, we have everything in growing order of coordinates, but on - strand it is the opposite.
    So we need to explain the tool how to deal with those fragments when extracting them.

    The "Region" tag is necessary to name each range differently, so that the _extract_one() helper can use these tags as key when
    building the dictionary of coverages per exonic fragment.
    """

    # First, convert object to a dataframe, for manipulation
    my_ranges = transcript_ranges.as_df()[['Chromosome', 'Feature', 'Start', 'End', 'transcript_id', 'gene_name', 'exon_number', 'exon_id', "Strand"]]

    # Check strandedness of the transcript features is unique (It should be, but errors in annotation happen).
    if my_ranges["Strand"].nunique() != 1:
        raise AssertionError("Transcript should be on one strand, but the selected one seems to have elements on both/none strands.")
    else:
        strandedness = my_ranges["Strand"].unique()

    # Extract features coordinates, divided between CDS, utr5 and utr3
    utr5_ranges_df = my_ranges[my_ranges["Feature"]== "five_prime_utr"].copy()
    cds_ranges_df = my_ranges[(my_ranges["Feature"]== "CDS") | (my_ranges["Feature"]== "stop_codon") ].copy()
    utr3_ranges_df = my_ranges[my_ranges["Feature"]== "three_prime_utr"  ].copy()

    # Resets the index, based on the order of features
    utr5_ranges_df = utr5_ranges_df.sort_values(by = "Start", ascending = True).reset_index(drop = True)
    cds_ranges_df = cds_ranges_df.sort_values(by = "Start", ascending = True).reset_index(drop = True)
    utr3_ranges_df = utr3_ranges_df.sort_values(by = "Start", ascending = True).reset_index(drop = True)

    # Applies the index to the "Region" tag (used by _extract_one() later)
    utr5_ranges_df["Region"] = utr5_ranges_df["Feature"] + "|" + utr5_ranges_df.index.astype(str)
    cds_ranges_df["Region"] = cds_ranges_df["Feature"] + "|" + cds_ranges_df.index.astype(str)
    utr3_ranges_df["Region"] = utr3_ranges_df["Feature"] + "|" + utr3_ranges_df.index.astype(str)

    # Output
    # If returning the coordinates as a single dataframe, it's important to order them based on strand they are on
    if return_individual_regions == False:

        if strandedness == "+":
            all_ranges = pd.concat([utr5_ranges_df, cds_ranges_df, utr3_ranges_df])
        elif strandedness == "-":
            all_ranges = pd.concat([utr3_ranges_df, cds_ranges_df, utr5_ranges_df,])

        assert all_ranges["Start"].is_monotonic_increasing, "Something when wrong when concatenating utrs and cds coordinates based on polarity of the strand. Check regions_dataframe_formatter"

        return all_ranges, strandedness

    # Returns coordinates as separate dataframes for utr5, cds, and utr3
    elif return_individual_regions == True:

        return utr5_ranges_df, cds_ranges_df, utr3_ranges_df, strandedness


# Function to derive lengths of CDS and UTRs from ranges
def extract_features_length(all_ranges):
    """

    :param all_ranges: A dataframe with the ranges for each region (CDS, utrs)
    :return: the length of the three regions
    """
    #Calculate length of each individual feature
    all_ranges["length"] = all_ranges["End"] - all_ranges["Start"]

    # Calculate total length of mRNA
    mRNA_length = all_ranges["length"].sum()

    # Calculate length by grouped feature (CDS, utrs)
    features_breakdown = all_ranges.groupby("Feature")["length"].sum()
    cds_length = features_breakdown["CDS"] + features_breakdown["stop_codon"]
    utr5_length = features_breakdown["five_prime_utr"]
    utr3_length = features_breakdown["three_prime_utr"]

    # Sanity check
    assert cds_length + utr5_length + utr3_length == mRNA_length

    return utr5_length, cds_length, utr3_length


# Function to extract coverage for a single transcript from a single bigwig
def extract_mRNA_regions_from_bigwig(all_ranges_df, bigwig_file):
    """
    :param all_ranges_df: dataframe with genomic ranges (coordinates) for all elements of the mRNA
    :param bigwig_file: the path to a BigWig file
    :return: full transcript coverage, as a dataframe. SRR is kept a column name

    This function uses the internal helper _extract_one to extract
    the coverage for each element (exon or exon fragment) making up the mRNA regions
    (utr5_tuple, cds_tuple, utr3_tuple).

    The result of the internal helper is parsed to create the full mRNA coverage dataframe.

    An assert statement makes sure the concatenation respected the sequence order, by checking
    the genomic coordinates associated with each position.

    """

    # Extract from a single bigwig the coverage for all regions whose coordinates are specified in the three dataframes
    all_regions_tuple = _extract_one(bigwig_file, all_ranges_df)

    # Concatenate the fragments into the full mature mRNA
    full_transcript = pd.concat(all_regions_tuple[1].values(), ignore_index=True)

    # Sanity check: if the dataframes have been correctly created, we expect the positions to be monotonic crescent
    assert full_transcript["Position"].is_monotonic_increasing, "Genomic coordinates are not in ascending order. Something went wrong"

    return full_transcript


# Function to extract the single transcript from multiple bigwigs
def extract_from_multiple_bigwigs(bigwig_files, all_ranges, strandedness, n_jobs = 1):
    """
    :param bigwig_files: List of paths pointing to bigwigs to use
    :param transcript_ranges: pyranges object from the .gtf annotation file
    :param n_jobs: how many cores to allocate
    :return: tuple. Study ID and dictionary of the three regions dataframes.
    """

    # Get the parent interpreter context to enable parallel execution
    # context  = mp.get_context("fork")

    # Extracts in parallel the coverages from each bigwig
    # NEEDS TO BE FIXED SO TO DIVIDE BIGWIGS INTO CHUNKS AND ASSIGN A CHUNK TO EACH WORKER!!
    results = Parallel(n_jobs=n_jobs,
                       backend="multiprocessing") (
        # Needs to be modified so to handle fused df...?
        delayed(extract_mRNA_regions_from_bigwig)(all_ranges, bw) for bw in bigwig_files
    )

    ### Need to define how to merge dataframes together, but it should be easy to do...
    all_studies_df = results[0]
    for df in results[1:]:
        all_studies_df = all_studies_df.merge(df, on="Position", how="outer")

    ### and flip here if the strand is - ... (we can flip all dfs at this point.
    if strandedness == "-":
        all_studies_df = all_studies_df.sort_values("Position", ascending=False).reset_index(drop=True)

    all_studies_df.fillna(0, inplace=True)

    return all_studies_df


# Introducing smoothing function:
def smoother(df, k=3):
    """
    Calculates the frame over the global aggregate, then smoothes by it.
    Also adds column "X" used by plotting function "subplot_single_transcript_plot"

    :param df: pandas.DataFrame
        DataFrame containing a 'Global Aggregate' column
    :param k: k-mer size to smooth the data. Default is 3.
            If the transcript is > 5000nt, increases smoothing window proportionally
    :return df: pandas.DataFrame
        DataFrame with the smoothed Global Aggregate
    """
    # Copies the position to column "X". Used to calculate Frame and for visualisation later
    df["X"] = df.index
    # Assign frame to each position
    df["Frame"] = (df["X"] % 3) + 1

    # Code to smooth cov by frame
    transcript_len = len(df.index)
    if transcript_len < 3500:
        k = 3
    elif transcript_len < 6000:
        k = 5
    elif transcript_len < 10000:
        k = 9
    else:
        k = 13

    k2 = k * 3
    for i in range(0, len(df), k2):
        # Select the slice of 9 rows using iloc
        chunk = df.iloc[i:i + k2]
        # out of those, groupby Frame, calculate the mean
        chunk_averages = chunk.groupby('Frame')['Global Aggregate'].transform("mean")
        # store the new values as a new column
        df.loc[chunk.index, "Global_Aggregate_smoothed"] = chunk_averages

    return df, k

# Helper. Adds a vertical line to a plot, with a label on top.
def _add_vline_and_label(p, x, max_y, text):
    """
    Helper function. Adds a vertical line and a label to the given Bokeh plot.
    Used to mark start/stop codons on single transcript plots

    :param p: bokeh.plotting.Figure
    :param x : int/float - x position for the line and label
    :param y_position : int/float - y position for the label
    :param text : str - label for the line

    """
    vline = Span(
        location=x,
        dimension="height",
        line_color="black",
        line_width=2,
        line_dash="dashed"
    )
    label = Label(
        x=x,
        y= max_y + max_y * 0.5,
        text=text,
        text_font_size="11pt",
        text_color="black",
        text_align="center",
        x_offset=0,
        y_offset= -20,
        y_units = "data",
        background_fill_color="white",
        background_fill_alpha=1.0
    )
    p.add_layout(vline)
    p.add_layout(label)


# New function to make Single transcript plot. Shows also frame info.
def subplot_single_transcript_plot(df, y_axis, x_range, k, start=0, stop=0, txid = ""):
    """
    Creates the single transcript plot, each frame coloured differently.

    :param df: -> input dataframe. Must have columns "Frame", "X" (position over transcript).
    :param y_axis: -> str, the name of the column with the data to be plotted on y axis
    :param x_range: -> x values, a factorRange object from Blokeh (to put different plots on same X together).
    :param start: -> int, position of START codon.
    :param stop: -> int, position of STOP codon.

    """
    # Ensure 'frame' is str for color mapping and legend
    df['Frame'] = df['Frame'].astype(str)

    # Color map to match R code
    frame_color_map = {
        "1": "#48cf4b",
        "2": "#68ace3",
        "3": "#e66343"
    }

    max_y = df[y_axis].max()

    p = figure(
        x_range=x_range,
        y_range = ( -1* df[y_axis].quantile(0.10), max_y * 1.5),
        sizing_mode="stretch_both",
        background_fill_color="white",
        x_axis_label="Position on the transcript (nt)",
        y_axis_label="Normalised Aggregate Counts",
        title= f"Global Aggregate for {txid}, smoothed by k = {k}"
    )

    # Remove gridlines and set white background (like theme_minimal + tweaks)
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.outline_line_color = None

    # Plot each frame as an area (overlay, not stacked)
    for frame, group in df.groupby('Frame'):
        group_sorted = group.sort_values('X')
        source = ColumnDataSource(group_sorted)
        p.varea(
            x="X", y1=0, y2=y_axis,
            source=source,
            fill_color=frame_color_map.get(frame, "#cccccc"),
            fill_alpha=0.6,
            legend_label=f"Frame {frame}",
        )

        p.vbar(
            x="X",
            top=y_axis,
            bottom=0,
            width=1,
            source=source,
            fill_color=frame_color_map.get(frame, "#cccccc"),
            fill_alpha=0,  # invisible until hover
            hover_fill_color=frame_color_map.get(frame, "#cccccc"),
            hover_fill_alpha=0.4,  # semi-transparent on hover
            line_color=None
        )

    # Style axes text (like axis.text.y = element_text(size=13))
    p.yaxis.axis_label_text_font_size = "13pt"
    p.xaxis.axis_label_text_font_size = "13pt"
    p.yaxis.major_label_text_font_size = "13pt"
    p.xaxis.major_label_text_font_size = "13pt"

    # Set legend title and placement
    p.legend.title = "Frame"
    p.legend.location = "top_center"
    p.legend.orientation = "horizontal"
    p.legend.click_policy = "hide"

    hover = HoverTool(
        tooltips=[
            ("Position", "@X"),
            ("Frame", "@Frame"),
            ("Value", f"@{y_axis}"),
        ],
        # Use 'vline' to snap tooltip to the nearest x-value (optional)
        mode='vline'
    )
    p.add_tools(hover)
    p.xaxis.major_label_text_font_size = '0pt'
    p.xaxis.major_tick_line_color = None
    p.xaxis.minor_tick_line_color = None

    # Define ranges to limit plot pan
    p.y_range.bounds = ( -1* df[y_axis].quantile(0.10), max_y*1.5)
    max_x = max(map(int, x_range.factors)) * 1.5
    min_x = min(map(int, x_range.factors)) - max_x * 0.2
    p.x_range.bounds = (min_x, max_x)

    # Add lines to highlight Start/Stop codons
    _add_vline_and_label(p, start, df[y_axis].max() ,"Start")
    _add_vline_and_label(p, stop, df[y_axis].max()  ,"Stop")

    return p

# Helper - Functions to create ORFplot
def get_sequence(fasta, df, strandedness):
    full_seq = ""

    for row in df.itertuples():
        chrom = row.Chromosome
        start = row.Start
        end = row.End
        seq = fasta[chrom][start:end].seq
        full_seq = full_seq + seq

    if strandedness == "+":
        return full_seq

    elif strandedness == "-":
        full_seq = str(Seq(full_seq).reverse_complement())
        return full_seq

# Helper - Map codons
def _map_type(codon):
    if codon == "ATG":
        return "ATG"
    elif codon in ["CTG", "GTG", "TTG"] : #, "GTG", "TTG"]:
        return "CTG"
    elif codon in ["TGA", "TAA", "TAG"]:
        return "stop codon"
    else:
        return None

# Helper - finds codon in a given frame
def _find_codons(seq, frame, codons):
    rows = []
    for i in range(frame - 1, len(seq) - 2, 3):
        codon = str(seq[i:i + 3])
        # if codon in codons:  # This way everything is mapped, regardless of given list
        rows.append({
            "codon": codon,
            "position": i,
            "Frame": int(frame),
            "type": _map_type(codon)
        })
    locations = pd.DataFrame(rows)

    return locations


def compute_codons_by_frame(sequence):
    codons = ["ATG", "CTG", "GTG", "TTG", "TGA", "TAA", "TAG"]
    frame_1 = _find_codons(sequence, 1, codons)
    frame_2 = _find_codons(sequence, 2, codons)
    frame_3 = _find_codons(sequence, 3, codons)

    temp = pd.concat([frame_1, frame_2, frame_3], ignore_index=True)
    all_reading_frames = temp.sort_values(by="position", ascending = True).reset_index(drop = True)

    return all_reading_frames

# Function to make one ORFplot for a given frame
def ORF_plotter(df, frame, x_range=None):
    # Masking as in your code
    plot_df = df.copy()

    mask = plot_df['Frame'] != frame
    plot_df.loc[mask, 'type'] = None

    # Convert to strings for categorical axes
    plot_df['type'] = plot_df['type'].astype(str)
    plot_df['position'] = plot_df['position'].astype(str)
    if x_range is None:
        x_range = plot_df['position'].tolist()

    frame_color_map = {
        1: "#48cf4b",
        2: "#68ace3",
        3: "#e66343"
    }

    bg_colour = frame_color_map[frame]

    # Bokeh requires a color list for the bars, so we map colors
    codon_color_map = {
        "ATG": "white",
        "stop codon": "black",
        "CTG": "#e3e3e3",
        "None": bg_colour,
        "nan": bg_colour,  # in case you get 'nan' from None/NaN conversion
    }
    # Fill bar colors accordingly
    plot_df['color'] = plot_df['type'].map(lambda x: codon_color_map.get(x, bg_colour))

    # Set up ColumnDataSource
    source = ColumnDataSource(plot_df)

    # Figure setup
    p = figure(
        x_range=x_range,
        y_range=(0, 1),
        height=75, width=700, title=f"Frame {frame}",
        toolbar_location="right",
        tools="wheel_zoom, reset"
    )
    bars = p.vbar(
        x='position',
        top=1,
        width=0.9,
        color='color',
        source=source,
        name="bars"
    )

    # Add interactive hover tooltip
    hover = HoverTool(
        tooltips=[
            ("Position", "@position"),
            ("Codon Type", "@type"),
            # Add other fields as you want, e.g. ("Value", "@value")
        ],
        renderers=[bars]
    )
    pan = PanTool(dimensions="width")
    p.add_tools(pan)
    p.toolbar.active_drag = pan
    p.add_tools(hover)

    p.y_range.bounds = (0, 1)

    # Hide y-axis
    p.yaxis.visible = False
    p.xaxis.visible = False
    p.yaxis.major_label_text_font_size = '0pt'
    p.yaxis.minor_tick_line_color = None
    p.yaxis.major_tick_line_color = None

    # Hide grid lines if desired
    p.ygrid.grid_line_color = None

    return p


# Make global aggregate plot (used in TranscriptID mode)
def plot_global_aggregate_with_zero(df, utr5_length = 0, cds_length = 0, txid = "", title="Global Aggregate of Counts for "):
    """
    Create an interactive plot with better control over axis ranges.

    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing a 'Global Aggregate' column
    utr5_length : int, optional
        Length of 5 UTR
    cds_length : int, optional
        Label of CDS
    title : str, optional
        Title for the plot

    Returns:
    --------
    plotly.graph_objects.Figure
        Interactive plotly figure with fixed axis behavior
    """

    if 'Global Aggregate' not in df.columns:
        raise ValueError("DataFrame must contain a 'Global Aggregate' column")

    fig = go.Figure()

    color_list = ["#d4ad07", "#21963b", " #b50000", "#3F4EC7", "#8737DA", "#1690D2"]
    chosen_color = random.choice(color_list)

    fig.add_trace(go.Scatter(
        x=df.index,
        y=df['Global Aggregate'],
        mode='lines',
        line=dict(width=2, color = chosen_color),
        hovertemplate='<b>Position:</b> %{x}<br><b>Total Counts:</b> %{y}<extra></extra>'
    ))

    fig.add_vline(
    x=utr5_length,  # your x-value
    line_dash="dash",  # dashed line
    line_color="white",  # or any color you want
    line_width=2
    )

    fig.add_annotation(
    x=utr5_length,
    y=1,  # or pick a y-value (if yref="paper", 1 is the top of the plot)
    yref="paper",
    text="CDS start",
    showarrow=False,
    font=dict(size=12, color="black"),
    bgcolor="white",
    bordercolor="black",
    borderpad=4
    )

    fig.add_vline(
    x=utr5_length+cds_length,  # your x-value
    line_dash="dash",  # dashed line
    line_color="white",  # or any color you want
    line_width=2
    )

    fig.add_annotation(
    x=utr5_length+cds_length,
    y=1,  # or pick a y-value (if yref="paper", 1 is the top of the plot)
    yref="paper",
    text="CDS end",
    showarrow=False,
    font=dict(size=12, color="black"),
    bgcolor="white",
    bordercolor="black",
    borderpad=4
    )

    fig.update_layout(
        title=dict(text=title + txid, x=0, font=dict(size=16)),
        xaxis_title="Position on transcript",
        yaxis_title="Total Counts",
        hovermode='x',
        showlegend=False,
        width=800,
        height=500,
        template='plotly_white',
        yaxis=dict(
            range=[0, df['Global Aggregate'].max() * 1.25],
            fixedrange=False,
            zeroline=True,
            zerolinecolor='rgba(0,0,0,0.5)',
            zerolinewidth=2
        ),
        xaxis=dict(
            rangeslider=dict(visible=False),
            type='linear',
            fixedrange=False
        )
    )

    return fig


def ranges_selector(full_mRNA_length, region):
    """
    A simple function to automatically render the selectors for ranges
    :param full_mRNA_length: Intuitively, the full cDNA length
    :param region: internal, needs to be different every time this function is called, or input keys will overlap and cause error
    :return: a tuple, with two int defining a range on the mRNA
    """
    x_start = st.number_input("Transcript region start coordinate (es: 10)",
                              min_value=0,
                              max_value=full_mRNA_length,
                              step=1,
                              key = f"start_{region}",
                              )
    x_end = st.number_input("Transcript region end coordinate (es: 100)",
                            min_value=0,
                            max_value=full_mRNA_length,
                            step=1,
                            key=f"end_{region}",
                            )
    st.write(f"Selected region {x_start} to {x_end}")

    return (x_start, x_end)