"""
home_helpers.py

This module contains functions used to render the Home page of the app, set up the analysis (organism) and format
inputs of the Upload mode. It works in sinergy with Home.py and relies on the functions in app_helpers.py

Author: Gionmattia Carancini
Date: 25-08-2025
"""

import os

import streamlit.elements.widgets.selectbox
import streamlit.runtime.session_manager

from app_helpers import *

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
            st.session_state.mappability_path = organism_config.get("mappability_path")

            bigwigs_dir = organism_config["bigwigs_dir"]
            st.session_state.fasta_path = organism_config["fasta_path"]
            st.session_state.bigwig_files = [
                os.path.abspath(os.path.join(bigwigs_dir, f))
                for f in os.listdir(bigwigs_dir)
                if os.path.isfile(os.path.join(bigwigs_dir, f))
            ]
            st.session_state.full_dataset_list = [
                f.replace(".bw", "") for f in os.listdir(bigwigs_dir)
                if os.path.isfile(os.path.join(bigwigs_dir, f))
            ]

            # These functions can be time-consuming
            st.session_state.fasta = read_fasta(st.session_state.fasta_path)
            st.session_state.gr = read_gtf_annotation(st.session_state.gtf_path)

        # 3. Set a flag indicating that the configuration is loaded
        st.session_state.config_loaded = True
        st.success(f"Loaded config for {st.session_state.selected_organism}!")


def genomic_coordinates_selection():
    """
    Function to render the UI for the "Genomic Coordinates" option in the Home page.
    """
    # Add a small explaination
    with st.expander("Show Genomic Mode Explaination"):
        st.markdown("""
        
        In this analysis mode, you are required to select three genomic intervals. 
        - The first one acts as your ***baseline***, to which you compare your... 
        - ...second one (usually a region you are interested in assessing its translational activity).
        - The third range should be a region closer to the second, but with ***no translational activity*** 
        (such as within a 3'trailer), that will be used as 'negative control' by the analysis.
        
        The resulting "RTA" (Ratio of Translational Activity) is a measure of how footprints densities is for the second range,
        compared to the baseline.
        
        ***NB1.*** As a general rule of thumb, try to select ranges of approximatively the same size and between 100bp - 2-3kb.
        
        ***NB2.*** This analysis mode is conceived to compare specific small genomic intervals. If you are interested
        in evaluating the translational activity for a transcript, please use the apposite mode.""")

    st.write("Please provide genomic ranges for the regions you want to compare: ")

    col1, _, _ = st.columns([3,2,1])
    with col1:
        genomic_range_1 = st.text_input("(Region 1 - Ex: *chr19:2,269,606-2,269,738*)")
        genomic_range_2 = st.text_input("(Region 2 - Ex: *chr19:2,271,786-2,271,949*)")
        genomic_range_n = st.text_input("(Noise Region - Ex: *chr19:2,273,160-2,273,406*)")
        # Save three ranges as a list
        genomic_ranges = [genomic_range_1, genomic_range_2, genomic_range_n]

    if st.button("Next"):

        # Generate ID for current run
        st.session_state.quirky_id = generate_quirky_id()
        # Initialize custom dataclasses (AnalysisStatus and GenomicCoordinatesMode) to keep track of analysis flow.
        st.session_state.analysis_runs[st.session_state.quirky_id] = AnalysisStatus(
            quirky_id=st.session_state.quirky_id,
            organism=st.session_state.selected_organism,
            mappability_path=st.session_state.mappability_path,
            genomic_mode=GenomicCoordinatesMode(
                active=True,
                genomic_ranges=genomic_ranges
            )
        )
        # Proceed to next page of the analysis
        st.switch_page("pages/metadata_filter.py")


def ensembl_txid_selection():
    """
    Function to render the UI for the "Ensembl Transcript ID" option in the Home page
    """
    st.write("Please provide a valid Ensembl stable transcript id  (*Ex. ENST00000602676, ENSMUST00000199678*): ")
    txid = st.text_input("Ensembl Transcript ID")

    if st.button("Next"):
        st.session_state.quirky_id = generate_quirky_id()
        # Create a new AnalysisStatus instance, associated to the quirky_id. Add flags for this mode.
        st.session_state.analysis_runs[st.session_state.quirky_id] = AnalysisStatus(
            quirky_id = st.session_state.quirky_id,
            organism = st.session_state.selected_organism,
            mappability_path= st.session_state.mappability_path,
            transcript_mode=TranscriptIDMode(
                active=True,
                txid = txid
            )
        )
        # Proceed to next page of the analysis
        st.switch_page("pages/metadata_filter.py")


def upload_selection():
    """
    Function to render the UI for the "File Upload" option in the Home page
    """

    st.write("Please select the type of input you want to load: ")
    if st.button("Genomic Coordinates"):
        st.session_state.tx_upload = False
        st.session_state.genomic_upload = True
    if st.button("Ensembl Transcript ID"):
        st.session_state.tx_upload = True
        st.session_state.genomic_upload = False

    # Block to work in "Genomic Upload Submode"
    if st.session_state.get("genomic_upload", False) is True:
        # Upload the three .csv files
        pre_feature_x_file  = st.file_uploader("Upload the .csv for the first region")
        post_feature_y_file = st.file_uploader("Upload the .csv for the second region")
        noise_n_file = st.file_uploader("Upload the .csv for the noise region")
        # Convert to dataframes
        try:
            pre_feature_x = pd.read_csv(pre_feature_x_file, index_col=0)
            post_feature_y = pd.read_csv(post_feature_y_file, index_col=0)
            noise_n = pd.read_csv(noise_n_file, index_col=0)
        except ValueError:
            st.stop()
        # Merge them into a single dataframe
        if all(dfs is not None for dfs in (pre_feature_x,
                                           post_feature_y,
                                           noise_n)):
            genomic_merged, genomic_merged_ranges = genomic_mode_formatter(pre_feature_x,
                                                                    post_feature_y,
                                                                    noise_n)

        if st.button("Next"):
            # Generate ID for current run
            st.session_state.quirky_id = generate_quirky_id()
            st.session_state.analysis_runs[st.session_state.quirky_id] = AnalysisStatus(
                quirky_id=st.session_state.quirky_id,
                organism=st.session_state.selected_organism,
                genomic_upload=st.session_state.genomic_upload,
                tx_upload=st.session_state.tx_upload,
                genomic_mode=GenomicCoordinatesMode(
                    active=True,
                    genomic_merged=genomic_merged,
                    genomic_merged_ranges=genomic_merged_ranges)
            )
            st.switch_page("pages/metadata_filter.py")

    # Block to work in "Transcript Upload Submode"
    if st.session_state.get("tx_upload", False) is True:
        uploaded_file = st.file_uploader("Upload the .csv with your transcript", type = "csv")
        try:
            all_studies_df = pd.read_csv(uploaded_file, index_col=0)
        except ValueError:
            st.stop()
        txid = st.text_input("Please name your Ensembl Transcript ID")
        # Allows to store the length of the three transcript regions
        utr5_length = st.number_input("Please, provide length of the 5' leader for selected transcript",
                        min_value=0,
                        step=1
                        )
        cds_length = st.number_input("Please, provide length of the main CDS for selected transcript",
                        min_value=0,
                        step=1
                        )
        utr3_length = st.number_input("Please, provide length of the 3' trailer for selected transcript",
                        min_value=0,
                        step=1
                        )
        tx_cov_extracted = True

        if st.button("Next"):
            # Generate ID for current run
            st.session_state.quirky_id = generate_quirky_id()
            st.session_state.analysis_runs[st.session_state.quirky_id] = AnalysisStatus(
                quirky_id=st.session_state.quirky_id,
                organism=st.session_state.selected_organism,
                genomic_upload=st.session_state.genomic_upload,
                tx_upload=st.session_state.tx_upload,
                transcript_mode=TranscriptIDMode(
                    active=True,
                    txid = txid,
                    all_studies_df=all_studies_df,
                    utr5_length=utr5_length,
                    cds_length=cds_length,
                    utr3_length=utr3_length,
                    tx_cov_extracted=tx_cov_extracted,
                    all_studies_df_backup = all_studies_df.copy(), # Generates back-up right away
                    all_studies_df_pre_aggregation_backup = all_studies_df.copy(), # Generates back-up right away
                )
            )
            st.switch_page("pages/metadata_filter.py")


##### Formatting Functions - These are input-validation functions #####
# Formats Genomic mode ranges
def genomic_mode_formatter(pre_feature_x:pd.DataFrame,
                           post_feature_y:pd.DataFrame,
                           noise_n:pd.DataFrame
                           ) -> tuple[pd.DataFrame, dict[str, tuple[int, int]]]:
    '''
    This module serves if the user selects upload - genomic. Basically, it formats the input so that they can be used
    by the same functions used by the "genomic coordinates" mode.

    Args:
        pre_feature_x: pd.DataFrame with coverages from genomic coordinates (rows are positions, columns are SRR Runs)
        post_feature_y: pd.DataFrame with coverages from genomic coordinates (rows are positions, columns are SRR Runs)
        noise_n: pd.DataFrame with coverages from genomic coordinates
        used as noise control (rows are positions, columns are SRR Runs)

    Returns:
        genomic_merged: pd.DataFrame, the three inputs, concatenated on top of each other.
        genomic_slices: dictionary, the coordinates to re-extract the original three DataFrames.
    '''

    # Create a single df for run_analysis()
    genomic_merged = pd.concat([pre_feature_x, post_feature_y, noise_n],
                               axis=0,
                               ignore_index=True)

    length_x = len(pre_feature_x)
    length_y = len(post_feature_y)
    length_n = len(noise_n)

    # Set up the dictionary of coordinates for run_analysis()
    genomic_slices = {"x_slice": (0, length_x),
                      "y_slice": (length_x,length_x + length_y),
                      "n_slice": (length_x + length_y, length_x + length_y + length_n)}

    return genomic_merged, genomic_slices

# Functions to validate the transcript upload should go down here, ideally.