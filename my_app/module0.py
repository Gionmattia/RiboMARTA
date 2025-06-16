"""
module0.py

This module contains all the functions used to extract the coverage values from a (genomic) bigwig.
It contains the functions to extract three uninterruped regions from the coordinates, as well as
the functions to extract a single transcript.
and to calculate the ratios.
Author: Gionmattia Carancini
Date: 26-05-2025
"""

# Import libraries
import pyBigWig
import pyranges as pr
import pandas as pd
from joblib import Parallel, delayed
import os
import multiprocessing as mp
mp.set_start_method("fork", force=True)
import plotly.graph_objects as go
import random
import re
import time

# # OLD VERSION - BACKUP Helper - Extracts the coverage info from a single bigwig
# def _extract_one(bw_path, regions_df):
#     """
#     :param bw_path: path to a signle .bigwig file
#     :param regions_df: dataframe with genomic regions to extract
#     :return: tuple. Study ID and dictionary of the three regions dataframes.
#     """
#
#     bw   = pyBigWig.open(bw_path)
#     sid  = os.path.basename(bw_path).split(".")[0]
#     out  = {}
#
#     # Convert the coordinates from 1-indexed inclusive to BigWig's 0-indexed end-exclusive [,)
#     # Not necessary, pyranges already converts the coordinates system
#     # regions_df["Start"] = regions_df["Start"] - 1
#
#     # iterate directly over the DataFrame
#     for row in regions_df.itertuples(index=False):
#         chrom, start, end, region = (
#             row.Chromosome, int(row.Start), int(row.End), row.Region
#         )
#         vals = bw.values(chrom, start, end, numpy=True)
#         out[region] = pd.DataFrame({
#             "Position": range(start, end),
#             sid: vals
#         })
#
#     bw.close()
#     return sid, out

# Helper - Extracts the coverage info from a single bigwig
def _extract_one(bw_path, regions_df):
    """
    :param bw_path: path to a signle .bigwig file
    :param regions_df: dataframe with genomic regions to extract
    :return: tuple. Study ID and dictionary of the three regions dataframes.
    """

    sid  = os.path.basename(bw_path).split(".")[0]
    out  = {}

    attempts = 0
    max_attempts = 2

    while attempts < max_attempts:
        try:
            bw = pyBigWig.open(bw_path)
            # iterate directly over the DataFrame
            for row in regions_df.itertuples(index=False):
                chrom, start, end, region = (
                    row.Chromosome, int(row.Start), int(row.End), row.Region
                )
                vals = bw.values(chrom, start, end, numpy=True)
                out[region] = pd.DataFrame({
                    "Position": range(start, end),
                    sid: vals
                })
            bw.close()
            break  # success, so exit the loop

        # Failsafe. Sometimes parallelisation cannot work if same file is opened at the same time by != forks.
        # We introduce a random sleep time for a worker, just to resolve the issue. If that does not work, then....
        # returns 0s instead.
        except:
            # Increase attempts counter
            attempts += 1
            if attempts < max_attempts:
                # Wait a random time to ensure no overlap
                wait = random.uniform(0.3, 1.5)
                time.sleep(wait)
            # If issue cannot be resolved by waiting, returns zeros.
            else:
                for row in regions_df.itertuples(index=False):
                    chrom, start, end, region = (
                        row.Chromosome, int(row.Start), int(row.End), row.Region
                    )
                    # Fill with zeros (or None if you prefer)
                    out[region] = pd.DataFrame({
                        "Position": range(start, end),
                        sid: None
                    })

    return sid, out


# Function to format the regions_df in a standardised way:

def format_input_df(my_list):
    """
    Unpacks genomic coordinates and formats them into a table that is used by the code to extract coverage

    :param my_list: A list of strings (genomic coordinates) as "chrx:56,741,761-56,761,390"
    :return: The dataframe with ranges to extract from genomic
    """
    data = []
    region_names = ["First Region", "Second Region", "Noise Region"]
    for i, reg in zip(my_list, region_names):
        fri = i.replace(",", "").replace(".", "").replace(":", "-").split("-")
        chrom = fri[0]
        if chrom.lower().startswith('chr'):
            chrom = chrom[3:]
        start = int(fri[1]) - 1  # Converts from 1-indexed closed to 0-index open
        end = int(fri[2])
        data.append({"Region": reg, "Chromosome": chrom, "Start": start, "End": end})
    return pd.DataFrame(data)
# Function to extract the coverage from all bigwigs and compile the three dataframes.

def extract_regions_from_bw_parallel(bigwig_files, regions_df, n_jobs=1):
    """
    This is the function used to extract coverage info when the user specifies the three genomic ranges manually.

    :param bigwig_files: list of paths to .bigwig files
    :param regions_df: dataframe with genomic regions to extract
    :param n_jobs: how many cores to allocate
    :return: three dataframes, with coverage from x, y, and n regions
    """

    # Extracts in parallel the coverages from each bigwig
    # SAME AS BELOW< NEEDS TO BE FIXED. BIGWIGS SHOULD BE GIVEN IN CHUNKS TO EACH WORKER
    results = Parallel(n_jobs=n_jobs,
                       backend="multiprocessing"
                       )(
        delayed(_extract_one)(bw, regions_df) for bw in bigwig_files
    )

    # Unpacks the results and compiles everything in a single dictionary
    final = {}
    for sample_id, bw_dict in results:
        for region, df_tmp in bw_dict.items():
            if region in final:
                final[region] = final[region].merge(
                    df_tmp, on="Position", how="outer"
                )
            else:
                final[region] = df_tmp

    pre_feature_x = final["First Region"]
    post_feature_y = final["Second Region"]
    noise_n = final["Noise Region"]

    return pre_feature_x, post_feature_y, noise_n

## From here are the functions to work with a single transcript.

# Function to load the ranges
def extract_tx_pyranges(gtf_path, txid):
    """
    Function to load the gtf annotation and

    :param gtf_path: path to the annotation file in use
    :param txid: the stable (ensebml) transcript id
    :return: a dataframe with ranges for each feature of txid
    """
    # Read the gtf file
    gr = pr.read_gtf(gtf_path)
    transcript_ranges = gr[gr.transcript_id == txid]
    if len(transcript_ranges) == 0:
        raise KeyError(f"transcript id {txid} was not found in the annotation file.")

    return transcript_ranges

# Function to format the transcript ranges, for extraction
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


# Function to extract coverage from bigwigs
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
def extract_from_multiple_bigwigs(bigwig_files, transcript_ranges, n_jobs = 1):
    """
    :param bigwig_files: List of paths pointing to bigwigs to use
    :param transcript_ranges: pyranges object from the .gtf annotation file
    :param n_jobs: how many cores to allocate
    :return: tuple. Study ID and dictionary of the three regions dataframes.
    """
    # Extract the ranges
    all_ranges, strandedness = regions_dataframe_formatter(transcript_ranges)

    # Get the parent interpreter context to enable parallel execution
    context  = mp.get_context("fork")

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


# Function to derive lengths of CDS and UTRs from ranges
def extract_features_length(all_ranges):
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


# Helpers
def global_aggregate(df, keep_pos = False):
    '''
    Creates a new dataframe with the "Global Aggregate" column, which is just the sum of counts per position
    '''

    df =  df.drop(columns="Position").copy()
    df.fillna(0, inplace=True)
    ga_df = df / df.sum()
    ga_df["Global Aggregate"] = ga_df.sum(axis = 1)
    if keep_pos == True:
        ga_df["Position"] = df["Position"]

    return ga_df

# Make plot
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