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
import plotly.express as px

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

# Organism Loader - Provides interface and code to load organism.
def organism_loader_interface():
    """
    Displays the Streamlit UI for selecting and loading an organism's data.

    This function creates a select box and a button. When the user clicks the
    button, it loads all necessary configuration and data files (metadata, paths, etc.)
    into st.session_state.

    Side Effects:
        - Populates st.session_state with data like 'gtf_path', 'metadata', 'fasta', etc.
        - Sets st.session_state.config_loaded = True upon successful loading.
    """


    # Load configuration file
    config = load_config("config.yaml")
    # Load path to annotation file, metadata, and the bigwigs in the app
    organisms_available = list(config.keys())

    # The selectbox determines which organism will be loaded
    selected_key = st.selectbox("Select from list", organisms_available)

    if st.button(f"Load {selected_key}"):
        # 1. Update session state with the selected organism and its specific config
        st.session_state.selected_organism = selected_key
        organism_config = config[st.session_state.selected_organism]

        # 2. Load all paths and data into session state
        with st.spinner(f"Loading data for {st.session_state.selected_organism}..."):
            st.session_state.gtf_path = organism_config["gtf_path"]
            st.session_state.metadata = pd.read_csv(organism_config["metadata_path"])
            bigwigs_dir = organism_config["bigwigs_dir"]
            st.session_state.fasta_path = organism_config["fasta_path"]
            st.session_state.bigwig_files = [
                os.path.abspath(os.path.join(bigwigs_dir, f))
                for f in os.listdir(bigwigs_dir)
                if os.path.isfile(os.path.join(bigwigs_dir, f))
            ]

            # These functions can be time-consuming
            st.session_state.fasta = read_fasta(st.session_state.fasta_path)
            st.session_state.gr = read_gtf_annotation(st.session_state.gtf_path)

        # 3. Set a flag indicating that the configuration is loaded
        st.session_state.config_loaded = True
        st.success(f"Loaded config for {st.session_state.selected_organism}!")

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
