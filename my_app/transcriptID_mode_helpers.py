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
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import Span, Label, FactorRange, HoverTool, PanTool

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
def regions_df_formatter(transcript_ranges):
    """
    :param transcript_ranges: a pyranges.pyranges_main.PyRanges object
    :param return_individual_regions: whether to return individual regions (3', 5', cds) dataframe coordinates or as a single one
    :return: tuple, a singe dataframe, with the coordinates for each exon fragment. The strandedness of the feature extracted.

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
    ranges_df = my_ranges[my_ranges["Feature"]== "exon"].copy()

    # Resets the index, based on the order of features
    ranges_df = ranges_df.sort_values(by = "Start", ascending = True).reset_index(drop = True)

    # Applies the index to the "Region" tag (used by _extract_one() later)
    ranges_df["Region"] = ranges_df["Feature"] + "|" + ranges_df["exon_number"]

    assert ranges_df["Start"].is_monotonic_increasing, "Something when wrong when extracting regions. Check regions_dataframe_formatter"

    return ranges_df, strandedness


# Function to derive lengths of CDS and UTRs from ranges
def extract_features_length(ranges, strandedness):
    """

    :param all_ranges: A dataframe with the ranges for each region (CDS, utrs)
    :return: the length of the three regions
    """

    # Extract Feature informations available.
    df = ranges.as_df()

    features = df["Feature"].unique()
    # Calculate length of each individual feature
    df["length"] = df["End"] - df["Start"]
    # Calculate length by grouped feature (CDS, utrs)
    features_breakdown = df.groupby("Feature")["length"].sum()

    # Retrieve lengths based on presence of features on annotation (defaults to 0 or None)
    cds_length   = features_breakdown.get("CDS", 0)
    stop_len  = features_breakdown.get("stop_codon", 0)
    utr5_len = features_breakdown.get("five_prime_utr", None)
    utr3_len = features_breakdown.get("three_prime_utr", None)

    # Combine stop codon to cds
    #cds_length = cds_len + stop_len

    # Fall back to coordinate-based UTR inference when annotation lacks explicit UTR features.
    # Handle each UTR independently — they can be missing separately.
    if "start_codon" in features and "stop_codon" in features:
        if utr5_len is None:
            utr5_len = (retrieve_5_pad(df, strandedness) if strandedness == "+"
                        else retrieve_3_pad(df, strandedness))
        if utr3_len is None:
            utr3_len = (retrieve_3_pad(df, strandedness) if strandedness == "+"
                        else retrieve_5_pad(df, strandedness))

    return utr5_len, cds_length, utr3_len


#  helper function for extract_features_lengths
def retrieve_5_pad(ranges_df, strandedness):
    df = ranges_df[~ranges_df["Feature"].isin(["transcript", "CDS", "five_prime_utr", "three_prime_utr"])].copy()

    if strandedness == "+":
        # Retrieve start position (target)
        target = df[df["Feature"] == "start_codon"]["Start"].iloc[0]
        # Extract all exons starting before the start codon
        relevant_exons = df[df['Start'] <= target].copy()
        # Use clip to have a "roof value"
        relevant_exons['Clipped_End'] = relevant_exons['End'].clip(upper=target)
        # retrieve position
        left_pad = (relevant_exons['Clipped_End'] - relevant_exons['Start']).sum()

    elif strandedness == "-":
        # Retrieve start position (target)
        target = df[df["Feature"] == "stop_codon"]["Start"].iloc[0]

        # Extract all exons starting before the start codon
        relevant_exons = df[df['Start'] <= target].copy()
        # Use clip to have a "roof value"
        relevant_exons['Clipped_End'] = relevant_exons['End'].clip(upper=target)
        # retrieve position
        left_pad = (relevant_exons['Clipped_End'] - relevant_exons['Start']).sum()

    return left_pad

#  helper function for extract_features_lengths
def retrieve_3_pad(ranges_df, strandedness):
    df = ranges_df[~ranges_df["Feature"].isin(["transcript", "CDS", "five_prime_utr", "three_prime_utr"])].copy()

    if strandedness == "+":
        # Retrieve start position (target)
        target = df[df["Feature"] == "stop_codon"]["End"].iloc[0]
        # Extract all exons starting before the start codon
        relevant_exons = df[df['End'] >= target].copy()
        # Use clip to have a "roof value"
        relevant_exons['Clipped_Start'] = relevant_exons['Start'].clip(lower=target)
        # retrieve position
        right_pad = (relevant_exons['End'] - relevant_exons['Clipped_Start']).sum()

    elif strandedness == "-":
        # Retrieve start position (target)
        target = df[df["Feature"] == "start_codon"]["End"].iloc[0]

        # Extract all exons starting before the start codon
        relevant_exons = df[df['End'] >= target].copy()
        # Use clip to have a "roof value"
        relevant_exons['Clipped_Start'] = relevant_exons['Start'].clip(lower=target)
        # retrieve position
        right_pad = (relevant_exons['End'] - relevant_exons['Clipped_Start']).sum()

    return right_pad


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


# Comment - Functionally speaking, does a bit too much... but not sure if I have time/will to unpack it.
def txid_data_extraction(current_run, bigwig_files, gr):
    """
    Function loads coverage information into a dataframe (all_studies_df) which is stored into current_run.
    Also extract features lengths (5', CDs, 3') and stores them in current_run object.

    :param current_run: current run - object where txid mode variables and results are stored. See app_helpers.py
    :param bigwig_files: absolute path to bigwig files.
    :param gr: annotation object, as loaded by pyranges.



    :return: all_ranges, strandedness
    """

    with st.spinner("I am retrieving the coverage data..."):

        # Internal Check. Stops execution if the txid is not found for any reason.
        try:
            # 1 Retrieve transcript ranges
            transcript_ranges = extract_tx_pyranges(gr, current_run.txid)
        except KeyError:
            st.error("I am sorry, but the Transcript ID given was not found in the annotation file."
                     "Please, restart from the Home page.")
            img_placeholder.image("martas/ko_marta.png", width=100)
            st.stop()

        # 2  Retrieve transcript features (CDS, UTRs) ranges and strandedness
        #all_ranges, strandedness = regions_dataframe_formatter(transcript_ranges)
        all_ranges, strandedness = regions_df_formatter(transcript_ranges)
        # 3  Retrieve features (CDS, UTRs) lengths
        #         current_run.utr5_length, current_run.cds_length, current_run.utr3_length = extract_features_length(
        #             all_ranges)
        current_run.utr5_length, current_run.cds_length, current_run.utr3_length = extract_features_length(
            transcript_ranges, strandedness)
        # 4 Extract coverage from given bigwigs
        current_run.all_studies_df = extract_from_multiple_bigwigs(bigwig_files,
                                                                   all_ranges,
                                                                   strandedness,
                                                                   n_jobs=1)
        # Updates logic flags within current_run
        current_run.loading_tx_cov = False
        current_run.tx_cov_extracted = True

    return all_ranges, strandedness


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


def data_aggregator(current_run, aggregation_option):

    current_run.all_studies_df = current_run.all_studies_df_pre_aggregation_backup.copy()

    if aggregation_option == None:
        current_run.meta_agg_option = None
        return

    elif aggregation_option  ==  "Aggregate by metadata field": # Metadata Aggregation

        # Consistency, remove position (though it won't be aggregated, but better safe than sorry)
        current_run.all_studies_df = current_run.all_studies_df.drop(columns = "Position")
        metadata = st.session_state.metadata.copy()
        mapping = metadata.set_index("Run")[current_run.meta_agg_option]
        current_run.all_studies_df = current_run.all_studies_df.T.groupby(mapping).sum().T
        current_run.all_studies_df["Position"] = current_run.all_studies_df_pre_aggregation_backup["Position"]

        return

    elif aggregation_option  == "Total Aggregation": # Total Aggregation

        current_run.all_studies_df = current_run.all_studies_df.drop(columns = "Position")
        current_run.all_studies_df = current_run.all_studies_df.sum(axis=1).to_frame("Global Aggregate")
        current_run.all_studies_df["Position"] = current_run.all_studies_df_pre_aggregation_backup["Position"]

        return

def bubble(text, bg_color="#000000", outline_color="#f0aa2d", text_color="#f0aa2d", outline_width=1):
    st.markdown(
        f'<div style="'
        f'background:{bg_color};'
        f'color:{text_color};'
        f'border:{outline_width}px solid {outline_color};'
        f'padding:10px 14px;'
        f'border-radius:10px;'
        f'display:inline-block;'
        f'font-weight:600;'
        f'">{text}</div>',
        unsafe_allow_html=True,
    )
