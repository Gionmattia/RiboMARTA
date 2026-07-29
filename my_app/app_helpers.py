"""
app_helpers.py

This module contains all the functions used by RiboMARTA in general, to manage states and workflows.
Author: Gionmattia Carancini
Date: 11-08-2025
"""

import random
import streamlit as st
import pandas as pd
from dataclasses import dataclass, field
import yaml
from pyfaidx import Fasta
import pyranges as pr
import plotly.graph_objects as go
from TranscriptID_results_helpers import multi_line_plot

# Words used to define the analysis ID. Need to be human-readable and quirky, so that the user remembers them.
ADJECTIVES = [
     "bicurious", "clever", "sleepy", "stressed", "caffeinated", "cozy", "furious", "quirky", "cheeky", "flirty", "questioning",
    "longing", "carnal", "communist", "turbocapitalist", "geeky",
    "accurate", "active", "adamant", "closeted", "curious", "sketchy",
    "calm", "eager", "fancy", "festive", "jolly", "kind", "lively", "silly", "zany", "brave",
]
NOUNS = [
    "Kant", "Plato", "Socrates", "Newton", "Michelangelo", "Galileo", "Spartacus", "Diogenes", "Caesar", "Augustus",
    "Descartes", "Pythagoras", "DaVinci", "Cleopatra", "Hannibal", "Pascal", "Loki", "Thor", "Achilles", "Odysseus"
]

# Function to generate the analysis ID. Will be stored in custom class.
def generate_quirky_id() -> str:
    adjective = random.choice(ADJECTIVES)
    noun = random.choice(NOUNS)
    quirky_id = f"{adjective}_{noun}"
    return quirky_id

# Import config file with paths to files used by RiboMARTA
@st.cache_resource
def load_config(path_to_yaml):
    with open(path_to_yaml, 'r') as f:
        return yaml.safe_load(f)

@st.cache_resource
def read_fasta(fasta_path):
    fasta = Fasta(fasta_path)
    return fasta

# Function to load the gtf file. Separated so to cache the gtf and speed-up re-executions.
@st.cache_data
def read_gtf_annotation(gtf_path):
    """
    This function loads a .gtf annotation file as a pyRanges object.
    Created so that data is cached and execution is faster on subsequent runs.

    :param gtf_path: path to annotation file
    :return: pyRanges object
    """
    gr = pr.read_gtf(gtf_path)
    #st.info(f".gtf took: {gtf_time:.2f} to load.")
    return gr


## Below we add the sub-classes used within AnalysisStatus
@dataclass
class MetadataFilters:
    """Dataclass to store the filtering criteria and selected bigwigs for the analysis"""
    filtering_criteria: dict = None
    selected_studies: list = None


@dataclass
class GenomicCoordinatesMode:
    """State for genomic coordinates analysis mode"""
    # Step 1 - Flag. Genomic Mode activated | User inputs stored.
    active: bool = False # same purpose as old "genomic_coordinates"
    genomic_ranges: list = None
    genomic_upload: bool = False

    # Step 2 - Genomic Ranges formatted as a dataframe.
    formatted_genomic_ranges: pd.DataFrame = None

    # Step 3 - Flag. Regions Extracted | Slices stored.
    regions_extracted: bool = False
    genomic_slices: list = None
    mappabilities: list = None
    genomic_merged: pd.DataFrame = None
    genomic_merged_backup: pd.DataFrame = None
    genomic_merged_ranges: dict = None
    genomic_merged_mappability_df: pd.DataFrame = None
    genomic_merged_mappability_ranges: dict = None
    smoothing_option: bool = False
    smoothing_criminal_thresh: int = 4
    aggregation_option: str = None  # Not yet implemented in genco!
    members_per_group : pd.DataFrame = None # Not yet implemented in genco!


    # Step 4 - Analysis Parameters stored.
    alpha_baseline: float = None
    alpha_target: float = None
    ci: int = None
    N: int = None
    n_jobs: int = None
    seed: int = None

    # Step 5 - Flag. Analysis done | Results and metadata stored.
    genco_analysis_done: bool = False
    results: pd.DataFrame = None
    results_meta: pd.DataFrame = None

    # Step 5b - Flag. BERTA Analysis done | Results and metadata stored.
    genco_BERTA_analysis_done: bool = False
    BERTA_results: pd.DataFrame = None
    BERTA_results_meta: pd.DataFrame = None


@dataclass
class TranscriptIDMode:
    # Step 1 - Flag. Txid mode activated | Ensembl Transcript ID.
    active: bool = False
    txid: str = ""

    # Step 2 - Flags. Coverage extraction and extracted -> used to manage flow/UI
    loading_tx_cov: bool = False
    tx_cov_extracted: bool = False

    # Step 3 - Annotation-derived metrics | Data Retrieved | Plots rendered
    utr5_length: int = None
    cds_length: int = None
    utr3_length: int = None
    all_studies_df: pd.DataFrame = None
    all_studies_df_backup: pd.DataFrame = None
    all_studies_df_pre_aggregation_backup: pd.DataFrame = None
    transcript_sequence: str = None
    mappability_df: pd.DataFrame = None
    tx_ga_plot: go.Figure = None
    triplet_ga_plot: go.Figure = None
    triplet_ga_plot_raw: go.Figure = None
    norm_triplet_ga_plot: go.Figure = None

    # Step 4 - User-selected regions | Analysis Parameters
    x_slice: tuple = None
    y_slice: tuple = None
    n_slice: tuple = None
    smoothing_option: bool = False
    aggregation_option: str = None
    members_per_group : pd.DataFrame = None
    meta_agg_option: str = None
    smoothing_criminal_thresh: int =  None
    txid_ranges: dict = None # Redundant?
    alpha_baseline: float = None
    alpha_target: float = None
    ci: int = None
    N: int = None
    n_jobs: int = None
    seed: int = None

    # Step 5 - Flag for legacy analysis | Legacy Results
    txid_analysis_done: bool = False
    txid_results: pd.DataFrame = None
    txid_results_meta: pd.DataFrame = None



# Dataclass to store the analysis, the logical flags, and the outputs. Each has its own quirky ID.
@dataclass
class AnalysisStatus:
    """
    Main analysis status container, reworked to use composition for readability.
    Each analysis has a unique quirky_id to encapsulate the run.

    """
    # General Info
    quirky_id: str = None
    organism: str = ""
    full_dataset_list: list = None
    mappability_path: str = None #Default is None

    # Studies selection
    filters: MetadataFilters = field(default_factory=MetadataFilters)

    # Analysis modes (only one active at a time)
    genomic_mode: GenomicCoordinatesMode = field(default_factory=GenomicCoordinatesMode)
    transcript_mode: TranscriptIDMode = field(default_factory=TranscriptIDMode)

    # Upload mode (they fall back into their respective modes)
    genomic_upload: bool = False
    tx_upload: bool = False

    # Shopping cart for SRR IDs
    flagged_datasets: list = field(default_factory=list)
    user_plots: list = field(default_factory=list)

    @property
    def active_mode(self) -> str:
        """Returns the analysis mode currently active"""
        if self.genomic_mode.active:
            return "genomic_coordinates"
        elif self.transcript_mode.active:
            return "ensembl_transcript_id"
        # elif self.upload_mode.active:
        #     return "upload"
        return "none"





# @dataclass
# class AnalysisStatus:
#     """
#     This custom class holds the analysis status. It comprises:
#     - Logical Flags
#     - User Inputs
#     - Outputs
#     Each analysis has a unique quirky_id, to encapsulate the analysis
#     """
#
#     # --- General Info: Analysis quirky ID and organism ---
#     quirky_id: str = None # to check with runs dictionary
#     organism: str = ""
#
#     # --- Parameters for the analysis (shared across modes) ---
#     ci: int = None
#     N: int = None
#     n_jobs: int = None
#
#     # --- Flags 1 - Analysis Mode ---
#     genomic_coordinates: bool = False
#     ensembl_transcript_id: bool = False
#     uploaded_file: bool = False
#
#     ## --- Variables and Flags for Genomic Coordinates Mode ---
#     # --- User Inputs for Genomic Mode ---
#     genomic_range_1: str = ""
#     genomic_range_2: str = ""
#     genomic_range_n: str = ""
#     genomic_ranges: list = None
#     formatted_genomic_ranges: pd.DataFrame = None
#     # --- Raw counts for Genomic Mode ---
#     genomic_slices: list = None
#     genomic_merged: pd.DataFrame = None
#     genomic_merged_ranges: dict = None
#     # --- Flags 2 - Extraction of regions ---
#     regions_extracted: bool = False
#
#     # --- Results Genomic Coordinates Mode ---
#     genco_results: pd.DataFrame = None
#     results_meta: pd.DataFrame = None
#     # --- Flag 3 - Analysis Ran for Genomic Coordinates mode ---
#     genco_analysis_done: bool = False
#     ### IMPORTANT: SO FAR, WE ARE NOT STORING RESULTS!! MAYBE THIS COULD BE CHANGED LATER ON...
#
#     ## --- Variables and Flags for Ensembl Transcript ID Mode ---
#     txid: str = ""
#
#     # --- Flags 3 - Tixd Mode ---
#     loading_tx_cov: bool = False
#     tx_cov_extracted: bool = False
#     show_ranges_selectors: bool = False
#
#     # --- Txid mode variables and parameters ---
#     utr5_length: int = None
#     cds_length: int = None
#     utr3_length: int = None
#     all_studies_df: pd.DataFrame = None
#     transcript_sequence: str = ""
#     tx_ga_plot: go.Figure = None
#     triplet_ga_plot: go.Figure = None
#     norm_triplet_ga_plot: go.Figure = None
#     # User-selected slices to compare in the analysis
#     x_slice: tuple = None
#     y_slice: tuple = None
#     n_slice: tuple = None
#     txid_ranges: dict = None
#
#     # --- Txid mode results ---
#     txid_results: pd.DataFrame = None
#     txid_results_meta: pd.DataFrame = None
#
#     # --- Flag 4 - Analysis Ran for Txid mode ---
#     txid_analysis_done: bool = False
#
#     # --- Variables and Flags for Update Mode ---
#     # Nb. This mode has not many flags. Theoretically, we want to converge the analysis into genomic or txid modes.
#     genomic_upload: bool = False
#     tx_upload: bool = False
#
#     # --- Flags for breakdowns ---#
#     tx_ranges = False
#     tx_regions = False
#
#     # --- Frame-aware Analysis Flags --- #
#     frame_aware: bool = False # If activated, will perform the RTA calculation while also comparing all possible Frames.
#     fa_txid_analysis_done: bool = False
#     # --- Frame-aware Analysis Inputs ---#
#     baseline_frame: int = None
#     fa_df: pd.DataFrame = None
#     # --- Frame-aware Analysis Outputs ---#
#     fa_txid_results: pd.DataFrame = None
#     fa_results_df: pd.DataFrame = None
#     fa_txid_results_meta: pd.DataFrame = None
#     fa_RTAs_plot = None
#
#     # --- Cart. A list where we store all SRR of interest ---#
    cart: list = field(default_factory=list)

# def compile_report(current_run):
#
#     # Get analysis mode
#     analysis_mode_flags = [current_run.genomic_coordinates,
#                            current_run.ensembl_transcript_id,
#                            current_run.uploaded_file]
#     flags_names = ["Genomic Coordinates", "Ensembl Transcript ID", "Uploaded File"]
#     analysis_mode = [name for name, value in zip(flags_names, analysis_mode_flags) if value]
#     analysis_mode = str(analysis_mode).replace("[", "").replace("]", "").replace("'","")
#
#     # Wraps in text - Still prints analysis_mode as a list... need to change it
#     text_opener = f"""
#     Analysis ID: {current_run.quirky_id} \n
#     Organism and reference: {current_run.organism} \n
#     Analysis Mode: {analysis_mode} \n
#     """
#
#     if analysis_mode == "Genomic Coordinates":
#         # Should I create a specific wrapper? like subfunctions?
#         user_input_genco_ranges = current_run.genomic_ranges
#         text_opener = text_opener + f"""Genomic Ranges selected: {user_input_genco_ranges}"""
#
#     elif analysis_mode == "Ensembl Transcript ID":
#         user_input_transcript_id = current_run.txid
#         ranges_selected = current_run.txid_ranges
#
#         txid_mode_info = f"""
#         Transcript selected: {user_input_transcript_id} \n
#         Baseline: {ranges_selected["x_slice"]} \n
#         Target: {ranges_selected["y_slice"]} \n
#         Noise region: {ranges_selected["n_slice"]} \n
#         """
#         text_opener = text_opener + txid_mode_info
#
#
#     return text_opener

## --- Sidebar elements test --- #
def srr_lookup_sidebar(data, current_run, analysis):

    st.sidebar.header("Quick Run (SRR) Lookup")
    srr_list = list(data.index)
    srr = st.sidebar.selectbox("Select Colour", srr_list)
    info = data.loc[srr]

    with st.sidebar.popover("Metadata", use_container_width=True):
        st.write(info)
    with st.sidebar.popover("Raw Data", use_container_width=True):
        # Add the plot (should just plot raw data for one selected SRR)
        srr_plot = fa_tx_line_plot(current_run.fa_df,
                                   current_run.utr5_length,
                                   current_run.cds_length,
                                   [srr],
                                   current_run.txid_ranges)

        st.plotly_chart(srr_plot, use_container_width=True)
    if st.sidebar.button("Add SRR to cart"):
        analysis.flagged_datasets.append(srr)
    if st.sidebar.button("Remove SRR from cart"):
        try:
            analysis.flagged_datasets.remove(srr)
        except:
            st.sidebar.write(f"{str(srr)} already NOT in your cart!")


## --- Sidebar elements for txid mode (could be fused to the one above, but later on) --- #
def srr_lookup_sidebar2(data, current_run, analysis):

    st.sidebar.header("Quick Run (SRR) Lookup")
    srr_list = list(data.index)
    srr = st.sidebar.selectbox("Select Colour", srr_list)
    info = current_run.txid_results_meta.loc[srr]

    with st.sidebar.popover("Metadata", use_container_width=True):
        st.write(info)

    with st.sidebar.popover("Raw Data", use_container_width=True):
        srr_plot = multi_line_plot(current_run.all_studies_df,
                                   current_run.utr5_length,
                                   current_run.cds_length,
                                   [srr],
                                   current_run.txid_ranges)

        st.plotly_chart(srr_plot, use_container_width=True)
    if st.sidebar.button("Add SRR to flagged ones"):
        if srr in analysis.flagged_datasets:
            st.sidebar.write(f"{str(srr)} has already been flagged!")
        else:
            analysis.flagged_datasets.append(srr)
            st.sidebar.write(f"{str(srr)} flagged!")
    if st.sidebar.button("Remove SRR "):
        try:
            analysis.flagged_datasets.remove(srr)
        except:
            st.sidebar.write(f"{str(srr)} NOT among flagged ones!")
    if analysis.flagged_datasets:
        st.sidebar.subheader("Flagged Datasets")
        for dataset in analysis.flagged_datasets:
            st.sidebar.write(dataset)
    else:
        st.sidebar.write("No datasets flagged! 🎉")



def quirky_id_checkpoint():
    if st.session_state.get("quirky_id", None) is None:
        col1, col2 = st.columns([1, 3])
        with col1:
            st.image("martas/error_marta.png", width=127)
        with col2:
            st.markdown("<br><br>", unsafe_allow_html=True)
            st.error("You haven't initialised any analysis yet, please return to Home page")
        st.stop()