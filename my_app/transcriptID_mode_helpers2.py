"""
transcriptID_mode_helpers2.py

This module contains functions for the transcriptID mode. Specifically, it contains functions used to render text, plots,
and other "front-end" options.
For functions relative to backend operations, such as extracting ranges or generating the plots in the first place,
please refer to transcriptID_mode_helpers.py

Author: Gionmattia Carancini
Date: 25-08-2025
"""
import streamlit as st
from app_helpers import AnalysisStatus
from transcriptID_mode_helpers import *


def quirky_id_checkpoint():
    if st.session_state.get("quirky_id", None) is None:
        col1, col2 = st.columns([1, 3])
        with col1:
            st.image("martas/error_marta.png", width=127)
        with col2:
            st.markdown("<br><br>", unsafe_allow_html=True)
            st.error("You haven't initialised any analysis yet, please return to Home page")
        st.stop()


def txid_mode_intro():
    """
    Renders Marta and a brief description of txid mode.
    :return:
    """
    col1, col2 = st.columns([4, 1])
    with col1:
        st.write("Welcome to the Transcript ID Analysis page. This analysis is broken down in two stages:")
        st.write("""1) By clicking on "Load Transcript", you will be shown a plot with the aggregate Ribo-Seq footprints 
        across all the available studies.""")
        st.write(
            "Since this first step can take time, you will be allowed to save the raw data if you wish so. Click on 'Dowload"
            "Data as CSV'")
        st.write(
            """2) Once the dynamic plot is visualised, you can click on 'Select Ranges' to specify the areas of the transcript
             you want to compare in the analysis.""")
    with col2:
        st.image("martas/explaining_marta.png", width=180)


def txid_start_checkpoint(current_run):
    """
    Checkpoint function. It checks if the user has submitted a transcriptID before allowing
    to further proceed in the analysis.
    Txid should be saved in the current_run (class: AnalysisStatus, see app_helpers.py)
    It overwrites the value of current_run.loading_tx_cov flag tp True (used for loading screen)
    :param current_run:
    :return:
    """

    txid = current_run.txid
    if txid == "":
        col1, col2 = st.columns([1, 3])
        with col1:
            st.image("martas/error_marta.png", width=127)
        with col2:
            st.markdown("<br><br>", unsafe_allow_html=True)
            st.error(" The 'Transcript ID' field was empty. Please provide a value before proceeding!")
    # Allows to extract the tx coverage
    else:
        if st.button("Load Transcript"):
            current_run.loading_tx_cov = True


def format_txid_df_download(current_run):
    """
    Sets up the df with coverage from all Runs for download.
    It adds the sequence information to it as well.
    :param current_run:
    :return:
    """
    # copy the current working dataframe
    export_dataframe = current_run.all_studies_df.copy()
    # Add the sequence to the copy, and allow for download.
    export_dataframe.insert(0, 'Sequence', list(current_run.transcript_sequence))
    # Allow to save the raw data as a .csv
    csv = export_dataframe.to_csv(index=True).encode('utf-8')
    return st.download_button(
        label="Download data as CSV",
        data=csv,
        file_name=f'{current_run.txid}_raw_counts.csv',
        mime='text/csv'
    )

