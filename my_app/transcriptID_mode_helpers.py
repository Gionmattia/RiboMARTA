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
import time # to time process execution

# Function to load the gtf file. Separated so to cache the gtf and speed-up re-executions.
@st.cache_data
def read_gtf_annotation(gtf_path):
    """
    This function loads a .gtf annotation file as a pyRanges object.
    Created so that data is cached and execution is faster on subsequent runs.

    :param gtf_path: path to annotation file
    :return: pyRanges object
    """
    start_gtf = time.time()
    gr = pr.read_gtf(gtf_path)
    end_t = time.time()
    gtf_time = end_t - start_gtf
    st.info(f".gtf took: {gtf_time:.2f} to load.")
    return gr

# Function to load the ranges
def extract_tx_pyranges(gtf_path, txid):
    """
    Function to load the gtf annotation and extract the ranges for selected transcript

    :param gtf_path: path to the annotation file in use
    :param txid: the stable (ensebml) transcript id
    :return: a dataframe with ranges for each feature of txid
    """
    # Read the gtf file
    gr = read_gtf_annotation(gtf_path)
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