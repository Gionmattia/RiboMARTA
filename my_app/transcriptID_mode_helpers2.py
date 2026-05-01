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
from mappability_helpers import bokeh_mappability_track
from mappability_helpers import txid_map_extraction
from bokeh.models import FactorRange
from bokeh.embed import file_html
from bokeh.resources import CDN
from bokeh.layouts import column
from streamlit_bokeh import streamlit_bokeh



def txid_mode_intro():
    """
    Renders Marta and a brief description of txid mode.
    :return:
    """
    col1, col2 = st.columns([4, 1])
    with col1:
        st.write("Welcome to the Transcript ID Analysis page.")
        st.write("""1) By clicking on "Load Transcript", you will be shown a plot with the aggregate Ribo-Seq footprints 
        across all the available studies. You are also allowed to download the data as a .csv, if you wish so.""")
        st.write(
            """2) Once the plots are visualised, you can click on 'Select Ranges' to specify the areas of the transcript
             you want to compare in the analysis. You will be asked to select a baseline, a target, and a noise region.""")
        st.write(
            """3) Finally, you can decide whether you want to take into consideration the region's reading frames 
            (frame-aware analysis) or not.""")
    with col2:
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.image("martas/explaining_marta.png", width=180)


def txid_start_checkpoint(current_run, analysis):
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
    # Allows to extract the tx coverage, if not in Upload Mode.
    else:
        if analysis.tx_upload == False:
            if st.button("Load Transcript"):
                current_run.loading_tx_cov = True
                # and also loads mappability track
                if analysis.mappability_path is not None:
                    mappability_df = txid_map_extraction(txid, analysis.mappability_path, st.session_state.gr)
                    mappability_df["X"] = mappability_df.index
                    other_col = [c for c in mappability_df.columns if c not in ['Position', 'X']][0]
                    current_run.mappability_df = mappability_df.rename(columns={other_col: 'Mappability'})


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


def ga_plots_UI_render(current_run, analysis):
    """
    Renders the plots stored (tx_ga_plot and triplet_ga_plot) as a st.expander.
    Also adds some info text.
    :param current_run: the class AnalysisStatus instance for the current analysis
    :return:
    """
    #with st.expander("Global Aggregate Plots"):
        # tab1, tab2, tab3 = st.tabs(["Simple Aggregate", "Raw Global Aggregate & ORF plot", "Normalised Global Aggregate & ORF plot"])
        # with tab1:
        #     st.write("#### Simple Global Aggregate plot")
        #     st.write("Below you can see a Single Transcript plot for your data.")
        #     st.write("The x-axis represents the position on the transcript, while the y-axis represents the Global Aggregate"
        #              "(the sum of raw counts from each study, in that position).")
        #     st.write("This plot is a simplified visualisation, with no information on ORFs or reading frames."
        #              "It's an easy-on-the eyes visualisation, less informative but useful for an initial exploration.")
            # if current_run.tx_ga_plot is not None:
            #     st.plotly_chart(current_run.tx_ga_plot, use_container_width=True)
    tab2, tab3 = st.tabs(["Raw Global Aggregate & ORF plot", "Normalised Global Aggregate & ORF plot"])
    with tab2:
        st.write("####  Raw Global Aggregate & ORF plot")
        st.write("Below you can see a Single Transcript plot for your data.")
        st.write("The x-axis represents the position on the transcript, while the y-axis represents the Global Aggregate"
                 "(the sum of raw counts from each study, in that position), smoothed by a window of size k for easy exploration.")
        st.write("Counts have been coloured by frame, and beneath the plot an ORF plot displaying starts and stops "
                 "has been added.")
        if current_run.triplet_ga_plot is not None:
            streamlit_bokeh(current_run.triplet_ga_plot, use_container_width=False)
            if st.button("Store in the report", key="triplet_ga_plot"):
                analysis.user_plots.append(current_run.triplet_ga_plot)
                st.success(f"Plot added! You currently have {len(analysis.user_plots)} plots.")

    with tab3:
        st.write("#### Normalised Global Aggregate & ORF plot")
        st.write("The x-axis represents the position on the transcript, while the y-axis represents the Normalised Global Aggregate"
                 "(each count in a position was normalised by the total counts of its study. The the normalised values have been summed to generate"
                 "the Normalised Global Aggregate), smoothed by a window of size k for easy exploration."
                 "The normalisation treats all studies as equivalent, regardless of the number of total counts."
                 "(ie. low-count studies contribute to the plot in the same way as high-count studies).")
        st.write("Counts have been coloured by frame, and beneath the plot an ORF plot displaying starts and stops "
                 "has been added.")
        if current_run.norm_triplet_ga_plot is not None:
            streamlit_bokeh(current_run.norm_triplet_ga_plot, use_container_width=False)


