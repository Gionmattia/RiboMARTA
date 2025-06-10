# collection of helpers to render home page of RiboMARTA
import streamlit as st
import pandas as pd
import yaml

# Import config file with paths to files used by RiboMARTA
@st.cache_resource
def load_config(path_to_yaml):
    with open(path_to_yaml, 'r') as f:
        return yaml.safe_load(f)

def genomic_coordinates_selection():
    """
    Function to render the UI for the "Genomic Coordinates" option in the Home page.
    ---
    State variables set: my_ranges
    State variables cleared: txid, picture, genomic_slices, all_studies_df, utr5_length,
        cds_length, utr3_length
    """
    st.write("Please provide genomic ranges for the regions you want to compare: ")

    col1, _, _ = st.columns([3,2,1])
    with col1:
        genomic_range_1 = st.text_input("(Region 1 - Ex: *chr1:75,440,360-75,440,569*)")
        genomic_range_2 = st.text_input("(Region 2 - Ex: *chr1:75,441,465-75,441,626*)")
        genomic_range_n = st.text_input("(Noise Region - Ex: *chr1:75,441,846-75,442,036*)")
        st.session_state.my_ranges = [genomic_range_1, genomic_range_2, genomic_range_n]

    if st.button("Next"):
        # session_state clean-up
        keys_to_reset = ["txid", "picture", "genomic_slices", "all_studies_df", "utr5_length",
        "cds_length", "utr3_length"]
        reset_session_state(keys_to_reset)
        # Proceed to next page of the analysis
        st.switch_page("pages/page_1.py")


def ensembl_txid_selection():
    """
    Function to render the UI for the "Ensembl Transcript ID" option in the Home page
    ---
    State variables set: txid
    State variables cleared: my_ranges, picture, genomic_slices, all_studies_df, utr5_length, cds_length, utr3_length
    """
    st.write("Please provide a valid Ensembl stable transcript id  (*Ex. ENST00000380152, ENSMUST00000199678*): ")
    st.session_state.txid = st.text_input("Ensembl Transcript ID")

    if st.button("Next"):
        # session_state clean-up
        keys_to_reset = ["my_ranges", "picture", "genomic_slices", "all_studies_df", "utr5_length",
         "cds_length", "utr3_length"]

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
        pre_feature_x = st.file_uploader("Upload the .csv for the first region")
        post_feature_y = st.file_uploader("Upload the .csv for the second region")
        noise_n = st.file_uploader("Upload the .csv for the noise region")
        # puts them into a single data-structure.
        st.session_state.genomic_slices = {"pre_feature": pre_feature_x, "post_feature": post_feature_y, "noise": noise_n}

    if st.session_state.tx_upload:
        # Allows to load the transcript dataframe
        uploaded_file = st.file_uploader("Upload the .csv with your transcript", type = "csv")
        try:
            st.session_state.all_studies_df = pd.read_csv(uploaded_file, index_col=0)
        except ValueError:
            st.stop()

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

        if st.session_state.tx_upload == True:
            keys_to_reset = ["my_ranges", "txid", "picture", "genomic_slices"]
            #Need to store transcript df and three lengths as a single unit
            reset_session_state(keys_to_reset)
        st.switch_page("pages/page_3.py")

# Helper - resets streamlit state variables
def reset_session_state(keys):
    """
    Clean-up function, helps remove session_state variables when switching between analysis options

    :param keys: list of session_state keys to remove
    :return:
    """
    for k in keys:
        st.session_state.pop(k, None)
