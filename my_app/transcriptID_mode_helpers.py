# transcriptID_mode_helpers.py
import pyfaidx
import pyranges
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


# Function to load the transcript_ranges
def extract_tx_pyranges(gr: pr.PyRanges, txid: str) -> pr.PyRanges:
    """
    Function to extract the transcript_ranges for a selected transcript out of the gr object.

    Parameters
    ----------
    gr: pr.PyRanges
        GTF annotation file, loaded as a PyRanges object.
    txid: str
        the stable (ensebml) transcript id.

    Returns
    -------
    transcript_ranges: pr.PyRanges
        a pyranges object with the genomic coordinates for any element of the selected transcript

    Raises
    ------
    KeyError
        if the txid is not found in the annotation file.
    """
    transcript_ranges = gr[gr.transcript_id == txid]
    if len(transcript_ranges) == 0:
        raise KeyError(f"transcript id {txid} was not found in the annotation file.")

    return transcript_ranges


# Function to format the transcript transcript_ranges, for extraction of coverage from BigWigs
def regions_df_formatter(transcript_ranges: pr.PyRanges) -> tuple[pd.DataFrame, str]:
    """
    Function to format the transcript_ranges for a selected transcript out of the gr object.
    The "Region" tag is necessary to name each range differently, so that the _extract_one() helper can use these tags as key when
    building the dictionary of coverages per exonic fragment.

    Parameters
    ----------
    transcript_ranges: pr.PyRanges
        a pyranges object with the genomic coordinates for any element of the selected transcript
    Returns
    -------
    tuple[pd.DataFrame, str]:
        The pd.DataFrame with the coordinates for each exon fragment of the given transcript.
        A string encoding the strand (+ or -) for the features extracted.
    Raises
    ------
    ValueError
        If the features to be extracted as a single transcript are annotated on no strand or both.
    """ 

    # First, convert object to a dataframe, for manipulation
    my_ranges = transcript_ranges.as_df()[['Chromosome', 'Feature', 'Start', 'End', 'transcript_id', 'gene_name', 'exon_number', 'exon_id', "Strand"]]

    # Check strandedness of the transcript features is unique (It should be, but errors in annotation happen).
    if my_ranges["Strand"].nunique() != 1:
        raise ValueError("Transcript should be on one strand, but the selected one seems to have elements on both/none strands.")
    else:
        strandedness = str(my_ranges["Strand"].iloc[0])

    # Extract features coordinates, divided between CDS, utr5 and utr3
    ranges_df = my_ranges[my_ranges["Feature"]== "exon"].copy()

    # Resets the index, based on the order of features
    ranges_df = ranges_df.sort_values(by = "Start", ascending = True).reset_index(drop = True)

    # Applies the index to the "Region" tag (used by _extract_one() later)
    ranges_df["Region"] = ranges_df["Feature"] + "|" + ranges_df["exon_number"]

    assert ranges_df["Start"].is_monotonic_increasing, "Something when wrong when extracting regions. Check regions_dataframe_formatter"

    return ranges_df, strandedness


# Function to derive lengths of CDS and UTRs from transcript_ranges
def extract_features_length(transcript_ranges: pr.PyRanges, strandedness: str) -> tuple[int,int,int]:
    """
    Function to derive the length of the features (CDS, 5' leader, 3' trailer) if annotated.
     If not, returns padding values (for visualisation purposes).

    Parameters
    ----------
    transcript_ranges: pr.PyRanges
        a pyranges object with the genomic coordinates for any element of the selected transcript
    strandedness: str
        A string encoding the strand (+ or -) for the features extracted.

    Returns
    -------
    utr5_len, cds_length, utr3_len: tuple[int, int, int]
        length of each region, if annotated.

    """
    # Extract Feature information available.
    df = transcript_ranges.as_df()

    features = df["Feature"].unique()
    # Calculate length of each individual feature
    df["length"] = df["End"] - df["Start"]
    # Calculate length by grouped feature (CDS, utrs)
    features_breakdown = df.groupby("Feature")["length"].sum()

    # Retrieve lengths based on presence of features on annotation (defaults to 0)
    cds_length   = features_breakdown.get("CDS", 0)
    stop_len  = features_breakdown.get("stop_codon", 0) # not needed
    utr5_len = features_breakdown.get("five_prime_utr", 0) # Was None
    utr3_len = features_breakdown.get("three_prime_utr", 0) # Was None

    # Backup: coordinate-based UTR when annotation lacks explicit UTR features.
    # Handle each UTR independently, as either of them can be missing separately.
    if "start_codon" in features and "stop_codon" in features:
        if utr5_len == 0:
            utr5_len = (retrieve_5_pad(df, strandedness) if strandedness == "+"
                        else retrieve_3_pad(df, strandedness))
        if utr3_len == 0:
            utr3_len = (retrieve_3_pad(df, strandedness) if strandedness == "+"
                        else retrieve_5_pad(df, strandedness))

    return utr5_len, cds_length, utr3_len


#  helper function for extract_features_lengths
def retrieve_5_pad(ranges_df: pd.DataFrame, strandedness:str) -> int:
    """
    Function to derive the padding for the 5' leader, based on strandedness

    Parameters
    ----------
    ranges_df: pd.DataFrame
        table with the genomic coordinates for any element of the selected transcript
    strandedness: str
        A string encoding the strand (+ or -) for the features extracted.

    Returns
    -------
    left_pad: int

    """
    df = ranges_df[~ranges_df["Feature"].isin(["transcript", "CDS", "five_prime_utr", "three_prime_utr"])].copy()

    if strandedness == "+":
        # Retrieve start position
        target = df[df["Feature"] == "start_codon"]["Start"].iloc[0]
        # Extract all exons starting before the start codon
        relevant_exons = df[df['Start'] <= target].copy()
        # Use clip to have a "roof value"
        relevant_exons['Clipped_End'] = relevant_exons['End'].clip(upper=target)
        # retrieve position
        left_pad = (relevant_exons['Clipped_End'] - relevant_exons['Start']).sum()

    elif strandedness == "-":
        # Retrieve start position
        target = df[df["Feature"] == "stop_codon"]["Start"].iloc[0]
        # Extract all exons starting before the start codon
        relevant_exons = df[df['Start'] <= target].copy()
        # Use clip to have a "roof value"
        relevant_exons['Clipped_End'] = relevant_exons['End'].clip(upper=target)
        # retrieve position
        left_pad = (relevant_exons['Clipped_End'] - relevant_exons['Start']).sum()

    return int(left_pad)


#  helper function for extract_features_lengths
def retrieve_3_pad(ranges_df: pd.DataFrame, strandedness: str) -> int:
    """
    Function to derive the padding for the 3' trailer, based on strandedness

    Parameters
    ----------
    ranges_df: pd.DataFrame
        table with the genomic coordinates for any element of the selected transcript
    strandedness: str
        A string encoding the strand (+ or -) for the features extracted.

    Returns
    -------
    right_pad: int

    """
    df = ranges_df[~ranges_df["Feature"].isin(["transcript", "CDS", "five_prime_utr", "three_prime_utr"])].copy()

    if strandedness == "+":
        # Retrieve start position
        target = df[df["Feature"] == "stop_codon"]["End"].iloc[0]
        # Extract all exons starting before the start codon
        relevant_exons = df[df['End'] >= target].copy()
        # Use clip to have a "roof value"
        relevant_exons['Clipped_Start'] = relevant_exons['Start'].clip(lower=target)
        # retrieve position
        right_pad = (relevant_exons['End'] - relevant_exons['Clipped_Start']).sum()

    elif strandedness == "-":
        # Retrieve start position
        target = df[df["Feature"] == "start_codon"]["End"].iloc[0]
        # Extract all exons starting before the start codon
        relevant_exons = df[df['End'] >= target].copy()
        # Use clip to have a "roof value"
        relevant_exons['Clipped_Start'] = relevant_exons['Start'].clip(lower=target)
        # retrieve position
        right_pad = (relevant_exons['End'] - relevant_exons['Clipped_Start']).sum()

    return int(right_pad)


# Function to extract coverage for a single transcript from a single bigwig
def extract_mRNA_regions_from_bigwig(all_ranges_df: pd.DataFrame, bigwig_file: str) -> pd.DataFrame:

    """
    This function uses the internal helper _extract_one to extract
    the coverage for each element (exon or exon fragment) making up the mRNA regions
    (utr5_tuple, cds_tuple, utr3_tuple).

    The result of the internal helper is parsed to create the full mRNA coverage dataframe.

    An assert statement makes sure the concatenation respected the sequence order, by checking
    the genomic coordinates associated with each position.

    Parameters
    ----------
    all_ranges_df: pd.DataFrame
        Data frame with the ranges for the transcript of interest.
    bigwig_file: str
        Path to a single bigwig file.

    Returns
    -------
    full_transcript: pd.DataFrame
        Data frame with the info for the transcript extracted from a single bigwig

    """
    # Extract from a single bigwig the coverage for all regions whose coordinates are specified in the three dataframes
    all_regions_tuple = _extract_one(bigwig_file, all_ranges_df)

    # Concatenate the fragments into the full mature mRNA
    full_transcript = pd.concat(all_regions_tuple[1].values(), ignore_index=True)

    # Sanity check: if the dataframes have been correctly created, we expect the positions to be monotonic crescent
    assert full_transcript["Position"].is_monotonic_increasing, "Genomic coordinates are not in ascending order. Something went wrong"

    return full_transcript


# Function to extract the single transcript from multiple bigwigs
def extract_from_multiple_bigwigs(bigwig_files: list,
                                  all_ranges: pd.DataFrame,
                                  strandedness: str,
                                  n_jobs:int = 1) -> pd.DataFrame:
    """
     This function extracts the coverage information from multiple bigwig files.

    Parameters
    ----------
    bigwig_files: list
        List of paths pointing to bigwigs to use, given as strings
    all_ranges: pd.DataFrame
        Data frame with the ranges for the transcript of interest.
    strandedness: str
        A string encoding the strand (+ or -) for the features extracted.
    n_jobs: int.
        how many cores to allocate to the extraction. Default is 1.

    Returns
    -------
    all_studies_df: pd.DataFrame
        Data frame with the raw count per position (rows) across all selected datasets (columns).

    """


    # Extracts in parallel the coverages from each bigwig
    # Nb. It would require to split the bigwig list across workers..which I did't do.
    # parallel implementation needs to be done later on!
    results = Parallel(n_jobs=n_jobs,
                       backend="multiprocessing") (
        # Needs to be modified so to handle fused df...?
        delayed(extract_mRNA_regions_from_bigwig)(all_ranges, bw) for bw in bigwig_files
    )

    ### Need to define how to merge dataframes together, but it should be easy to do...
    all_studies_df = results[0]
    for df in results[1:]:
        all_studies_df = all_studies_df.merge(df, on="Position", how="outer")

    ### and flip here if the strand is - ... (we can flip all dfs at this point).
    if strandedness == "-":
        all_studies_df = all_studies_df.sort_values("Position", ascending=False).reset_index(drop=True)

    all_studies_df.fillna(0, inplace=True)

    return all_studies_df


# def txid_data_extraction(current_run, bigwig_files, gr):
#
#     """
#     Wrapper function loads coverage information into a dataframe (all_studies_df) which is stored into current_run.
#     Also extract features lengths (5', CDs, 3') and stores them in current_run object.
#
#     Parameters
#     ----------
#     current_run: TranscriptIDMode
#         Custom DataClass instance storing all the info for the current analysis.
#     bigwig_files: list
#         List of paths pointing to bigwigs to use, given as strings
#     gr: pr.PyRanges
#         GTF annotation file, loaded as a PyRanges object.
#
#     Returns
#     -------
#     all_studies_df: pd.DataFrame
#         Data frame with the raw count per position (rows) across all selected datasets (columns).
#     strandedness: str
#         A string encoding the strand (+ or -) for the features extracted.
#     """
#
#     with st.spinner("I am retrieving the coverage data..."):
#
#         # Internal Check. Stops execution if the txid is not found for any reason.
#         try:
#             # 1 Retrieve transcript transcript_ranges
#             transcript_ranges = extract_tx_pyranges(gr, current_run.txid)
#         except KeyError:
#             st.error("I am sorry, but the Transcript ID given was not found in the annotation file."
#                      "Please, restart from the Home page.")
#             st.stop()
#
#         # 2  Retrieve transcript features (CDS, UTRs) transcript_ranges and strandedness
#         all_ranges, strandedness = regions_df_formatter(transcript_ranges)
#         # 3  Retrieve features (CDS, UTRs) lengths
#         current_run.utr5_length, current_run.cds_length, current_run.utr3_length = extract_features_length(
#             transcript_ranges, strandedness)
#         # 4 Extract coverage from given bigwigs
#         current_run.all_studies_df = extract_from_multiple_bigwigs(bigwig_files,
#                                                                    all_ranges,
#                                                                    strandedness,
#                                                                    n_jobs=1)
#         # Updates logic flags within current_run
#         current_run.loading_tx_cov = False
#         current_run.tx_cov_extracted = True
#
#     return all_ranges, strandedness


# Helper - Functions to create ORFplot
def get_sequence(fasta: pyfaidx.Fasta, df: pd.DataFrame, strandedness: str) -> str:
    """
     This function extracts the sequence information for a transcript from the genome FASTA file.

    Parameters
    ----------
    fasta: pyfaidx.Fasta
        Fasta file containing the genome for the selected organism
    df: pd.DataFrame
        Data frame with the ranges for the transcript of interest
    strandedness: str
        A string encoding the strand (+ or -) for the features extracted

    Returns
    -------
    full_seq: str
        The sequence for the transcript of interest

    """
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


def ranges_selector(full_mRNA_length: int, region: str) -> tuple[int,int]:
    """
    A simple renderer for numeric inputs. Used to allow user to select ranges for the regions of interest to test
    (baseline, target, noise control)

    Parameters
    ----------
    full_mRNA_length: int
        length of the full mRNA
    region: str
        the region of interest

    Returns
    -------
    (x_start, x_end): tuple(int,int)
        The ranges for the region of interest, as a tuple
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

# Should move elsewhere... It doesn't feel like belonging here.
def data_aggregator(current_run, aggregation_option):

    current_run.all_studies_df = current_run.all_studies_df_pre_aggregation_backup.copy()

    if aggregation_option == None:
        current_run.meta_agg_option = None
        return

    elif aggregation_option  ==  "Aggregate by metadata field": # Metadata Aggregation

        # Consistency, remove position (though it won't be aggregated, but better safe than sorry)
        current_run.all_studies_df = current_run.all_studies_df.drop(columns = "Position")
        metadata = st.session_state.metadata.copy()
        #subset to keep only selected meta
        sel_studies = current_run.all_studies_df.columns
        metadata = metadata[metadata["Run"].isin(sel_studies)]

        mapping = metadata.set_index("Run")[current_run.meta_agg_option]

        current_run.members_per_group = mapping.value_counts().to_frame()
        current_run.members_per_group["Group"] = current_run.members_per_group.index
        current_run.all_studies_df = current_run.all_studies_df.T.groupby(mapping).sum().T
        current_run.all_studies_df["Position"] = current_run.all_studies_df_pre_aggregation_backup["Position"]

        return

    elif aggregation_option  == "Total Aggregation": # Total Aggregation

        current_run.all_studies_df = current_run.all_studies_df.drop(columns = "Position")
        current_run.all_studies_df = current_run.all_studies_df.sum(axis=1).to_frame("Global Aggregate")
        current_run.all_studies_df["Position"] = current_run.all_studies_df_pre_aggregation_backup["Position"]

        return

