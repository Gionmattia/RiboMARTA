"""
app_helpers.py

This module contains all the functions used by RiboMARTA in general, to manage states and workflows.
Author: Gionmattia Carancini
Date: 11-08-2025
"""

import random
import streamlit as st
import pandas as pd
import os
from dataclasses import dataclass, field
from typing import List, Dict, Any
import yaml
from pyfaidx import Fasta
import pyranges as pr
import plotly.graph_objects as go
from frame_aware_analysis_helpers import fa_tx_line_plot
#from TranscriptID_results_helpers import multi_line_plot

# Words used to define the analysis ID. Need to be human-readable and quirky, so that the user remembers them.
ADJECTIVES = [
     "bicurious", "clever", "sleepy", "stressed", "caffeinated", "cozy", "furious", "quirky", "cheeky", "flirty", "questioning",
    "longing", "carnal", "communist", "turbocapitalist", "geeky"
    # "accurate", "active", "adamant", "closeted", "curious",
    # "calm", "eager", "fancy", "festive", "jolly", "kind", "lively", "silly", "zany", "brave",
]
NOUNS = [
    "Kant", "Plato", "Socrates", "Newton", "Michelangelo", "Galileo", "Spartacus",
    "Descartes", "Pythagoras", "DaVinci", "Cleopatra", "Hannibal", "Pascal", "Loki", "Thor"
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

# Dataclass to store the analysis, the logical flags, and the outputs. Each has its own quirky ID.
@dataclass
class AnalysisStatus:
    """
    This custom class holds the analysis status. It comprises:
    - Logical Flags
    - User Inputs
    - Outputs
    Each analysis has a unique quirky_id, to encapsulate the analysis
    """

    # --- General Info: Analysis quirky ID and organism ---
    quirky_id: str = None # to check with runs dictionary
    organism: str = ""

    # --- Parameters for the analysis (shared across modes) ---
    ci: int = None
    N: int = None
    n_jobs: int = None

    # --- Flags 1 - Analysis Mode ---
    genomic_coordinates: bool = False
    ensembl_transcript_id: bool = False
    uploaded_file: bool = False

    ## --- Variables and Flags for Genomic Coordinates Mode ---
    # --- User Inputs for Genomic Mode ---
    genomic_range_1: str = ""
    genomic_range_2: str = ""
    genomic_range_n: str = ""
    genomic_ranges: list = None
    formatted_genomic_ranges: pd.DataFrame = None
    # --- Raw counts for Genomic Mode ---
    genomic_slices: list = None
    genomic_merged: pd.DataFrame = None
    genomic_merged_ranges: dict = None
    # --- Flags 2 - Extraction of regions ---
    regions_extracted: bool = False

    # --- Results Genomic Coordinates Mode ---
    genco_results: pd.DataFrame = None
    results_meta: pd.DataFrame = None
    # --- Flag 3 - Analysis Ran for Genomic Coordinates mode ---
    genco_analysis_done: bool = False
    ### IMPORTANT: SO FAR, WE ARE NOT STORING RESULTS!! MAYBE THIS COULD BE CHANGED LATER ON...

    ## --- Variables and Flags for Ensembl Transcript ID Mode ---
    txid: str = ""

    # --- Flags 3 - Tixd Mode ---
    loading_tx_cov: bool = False
    tx_cov_extracted: bool = False
    show_ranges_selectors: bool = False

    # --- Txid mode variables and parameters ---
    utr5_length: int = None
    cds_length: int = None
    utr3_length: int = None
    all_studies_df: pd.DataFrame = None
    transcript_sequence: str = ""
    tx_ga_plot: go.Figure = None
    triplet_ga_plot: go.Figure = None
    # User-selected slices to compare in the analysis
    x_slice: tuple = None
    y_slice: tuple = None
    n_slice: tuple = None
    txid_ranges: dict = None

    # --- Txid mode results ---
    txid_results: pd.DataFrame = None
    txid_results_meta: pd.DataFrame = None

    # --- Flag 4 - Analysis Ran for Txid mode ---
    txid_analysis_done: bool = False

    # --- Variables and Flags for Update Mode ---
    # Nb. This mode has not many flags. Theoretically, we want to converge the analysis into genomic or txid modes.
    genomic_upload: bool = False
    tx_upload: bool = False

    # --- Flags for breakdowns ---#
    tx_ranges = False
    tx_regions = False

    # --- Frame-aware Analysis Flags --- #
    frame_aware: bool = False # If activated, will perform the RTA calculation while also comparing all possible Frames.
    fa_txid_analysis_done: bool = False
    # --- Frame-aware Analysis Inputs ---#
    baseline_frame: int = None
    fa_df: pd.DataFrame = None
    # --- Frame-aware Analysis Outputs ---#
    fa_txid_results: pd.DataFrame = None
    fa_results_df: pd.DataFrame = None
    fa_txid_results_meta: pd.DataFrame = None
    fa_RTAs_plot = None

    # --- Cart. A list where we store all SRR of interest ---#
    cart: list = field(default_factory=list)

def compile_report(current_run):

    # Get analysis mode
    analysis_mode_flags = [current_run.genomic_coordinates,
                           current_run.ensembl_transcript_id,
                           current_run.uploaded_file]
    flags_names = ["Genomic Coordinates", "Ensembl Transcript ID", "Uploaded File"]
    analysis_mode = [name for name, value in zip(flags_names, analysis_mode_flags) if value]
    analysis_mode = str(analysis_mode).replace("[", "").replace("]", "").replace("'","")

    # Wraps in text - Still prints analysis_mode as a list... need to change it
    text_opener = f"""
    Analysis ID: {current_run.quirky_id} \n
    Organism and reference: {current_run.organism} \n
    Analysis Mode: {analysis_mode} \n
    """

    if analysis_mode == "Genomic Coordinates":
        # Should I create a specific wrapper? like subfunctions?
        user_input_genco_ranges = current_run.genomic_ranges
        text_opener = text_opener + f"""Genomic Ranges selected: {user_input_genco_ranges}"""

    elif analysis_mode == "Ensembl Transcript ID":
        user_input_transcript_id = current_run.txid
        ranges_selected = current_run.txid_ranges

        txid_mode_info = f"""
        Transcript selected: {user_input_transcript_id} \n 
        Baseline: {ranges_selected["x_slice"]} \n
        Target: {ranges_selected["y_slice"]} \n
        Noise region: {ranges_selected["n_slice"]} \n
        """
        text_opener = text_opener + txid_mode_info


    return text_opener

## --- Sidebar elements test --- #
def srr_lookup_sidebar(data, current_run):

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
        current_run.cart.append(srr)
    if st.sidebar.button("Remove SRR from cart"):
        try:
            current_run.cart.remove(srr)
        except:
            st.sidebar.write(f"{str(srr)} already NOT in your cart!")
