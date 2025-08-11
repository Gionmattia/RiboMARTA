# collection of helpers to render home page of RiboMARTA
import streamlit as st
import pandas as pd
from module0 import genomic_mode_formatter
from app_helpers import generate_quirky_id
from app_helpers import AnalysisStatus

def genomic_coordinates_selection():
    """
    Function to render the UI for the "Genomic Coordinates" option in the Home page.
    ---
    State variables set: my_ranges
    State variables cleared: txid, picture, genomic_slices, all_studies_df, utr5_length,
        cds_length, utr3_length
    """
    # Add a small explaination
    with st.expander("Show Genomic Mode Explaination"):
        st.write("""
        
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
        #### THIS IS WHERE WE IMPLEMENT THE NEW ANALYSIS FLOW CONTROL...
        # Generate ID for current run
        st.session_state.quirky_id = generate_quirky_id()
        # Create a new Analysis Status class, associated to the quirky_id. Add flags for this mode.
        st.session_state.runs[st.session_state.quirky_id] = AnalysisStatus(quirky_id = st.session_state.quirky_id,
                                                            organism = st.session_state.selected_organism,
                                                            genomic_coordinates = True,
                                                            genomic_range_1 = genomic_range_1,
                                                            genomic_range_2=genomic_range_2,
                                                            genomic_range_n=genomic_range_n,
                                                            genomic_ranges = genomic_ranges)

        # session_state clean-up
        keys_to_reset = ["txid", "picture", "genomic_slices", "all_studies_df", "utr5_length",
        "cds_length", "utr3_length"]
        reset_session_state(keys_to_reset)
        # Reset regions_extracted (used to display page_1 breakdowns)
        st.session_state.regions_extracted = False
        # Proceed to next page of the analysis
        st.switch_page("pages/page_1.py")


def ensembl_txid_selection():
    """
    Function to render the UI for the "Ensembl Transcript ID" option in the Home page
    ---
    State variables set: txid
    State variables cleared: my_ranges, picture, genomic_slices, all_studies_df, utr5_length, cds_length, utr3_length
    """
    st.write("Please provide a valid Ensembl stable transcript id  (*Ex. ENST00000602676, ENSMUST00000199678*): ")
    txid = st.text_input("Ensembl Transcript ID")

    if st.button("Next"):
        st.session_state.quirky_id = generate_quirky_id()
        # Create a new Analysis Status class, associated to the quirky_id. Add flags for this mode.
        st.session_state.runs[st.session_state.quirky_id] = AnalysisStatus(quirky_id = st.session_state.quirky_id,
                                                            organism = st.session_state.selected_organism,
                                                            ensembl_transcript_id = True,
                                                            txid = txid)
        # session_state clean-up
        keys_to_reset = ["my_ranges", "picture", "genomic_slices","genomic_merged", "genomic_ranges", "all_studies_df", "utr5_length",
         "cds_length", "utr3_length", "tx_upload", "genomic_upload", "cov_extracted"]
        reset_session_state(keys_to_reset)
        # Proceed to next page of the analysis
        st.switch_page("pages/page_2.py")

def upload_selection():
    """
    Function to render the UI for the "File Upload" option in the Home page
    ---
    State variables set (for genomic_upload): genomic_slices
    State variables set (for tx_upload): all_studies_df, utr5_length, cds_length, utr3_length

    State variables cleared (for genomic_upload): my_ranges, txid, picture, all_studies_df,
                                                    utr5_length, cds_length, utr3_length
    State variables cleared (for tx_upload): my_ranges, txid, picture, genomic_slices
    """
    # Initialise the state variables
    if "genomic_upload" not in st.session_state:
        st.session_state.genomic_upload = False
    if "tx_upload" not in st.session_state:
        st.session_state.tx_upload = False

    st.write("Please select the type of input you want to load: ")

    if st.button("Genomic Coordinates"):
        st.session_state.genomic_upload = True
        st.session_state.tx_upload = False
    if st.button("Ensembl Transcript ID"):
        st.session_state.genomic_upload = False
        st.session_state.tx_upload = True

    if st.session_state.genomic_upload:
        # Allows to load three dataframes
        pre_feature_x_file  = st.file_uploader("Upload the .csv for the first region")
        post_feature_y_file = st.file_uploader("Upload the .csv for the second region")
        noise_n_file = st.file_uploader("Upload the .csv for the noise region")
        # Convert to dataframes
        pre_feature_x = pd.read_csv(pre_feature_x_file, index_col=0)
        post_feature_y = pd.read_csv(post_feature_y_file, index_col=0)
        noise_n = pd.read_csv(noise_n_file, index_col=0)
        # Create a single df and coordinates, to be used by run_analysis()


        # Patchy, but it works for now
        if all(dfs is not None for dfs in (pre_feature_x, post_feature_y, noise_n)):
            st.session_state.genomic_merged, st.session_state.genomic_ranges = genomic_mode_formatter(pre_feature_x,
                                                                                                  post_feature_y,
                                                                                                  noise_n)

    if st.session_state.tx_upload:
        # Allows to load the transcript dataframe
        uploaded_file = st.file_uploader("Upload the .csv with your transcript", type = "csv")
        try:
            st.session_state.all_studies_df = pd.read_csv(uploaded_file, index_col=0)
        except ValueError:
            st.stop()

        st.session_state.txid = st.text_input("Please name your Ensembl Transcript ID")

        # Allows to store the length of the three transcript regions
        st.session_state.utr5_length = st.number_input("Please, provide length of the 5' leader for selected transcript",
                        min_value=0,
                        step=1
                        )
        st.session_state.cds_length = st.number_input("Please, provide length of the main CDS for selected transcript",
                        min_value=0,
                        step=1
                        )
        st.session_state.utr3_length = st.number_input("Please, provide length of the 3' trailer for selected transcript",
                        min_value=0,
                        step=1
                        )

    if st.button("Next"):
        # session_state clean-up
        if st.session_state.genomic_upload == True:
            keys_to_reset = ["my_ranges", "txid", "picture" , "all_studies_df", "utr5_length",
                             "cds_length", "utr3_length"]
            reset_session_state(keys_to_reset)

            st.switch_page("pages/page_4.py")

        if st.session_state.tx_upload == True:
            keys_to_reset = ["my_ranges", "picture", "genomic_merged", "genomic_ranges"]
            #Need to store transcript df and three lengths as a single unit
            reset_session_state(keys_to_reset)
            st.session_state.cov_extracted = True
            st.switch_page("pages/page_2.py")

# Helper - resets streamlit state variables
def reset_session_state(keys):
    """
    Clean-up function, helps remove session_state variables when switching between analysis options

    :param keys: list of session_state keys to remove
    :return:
    """
    for k in keys:
        st.session_state.pop(k, None)
