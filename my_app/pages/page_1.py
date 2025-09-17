# page_1.py
from module0 import *
import streamlit as st
import pandas as pd
from home_helpers import genomic_mode_formatter

# We are retrieving the AnalysysStatus class corresponding to the quirky_id, from the dictionary "runs"
current_run = st.session_state.runs[st.session_state.quirky_id]
# This bit tests your current id matches the AnalysisStatus ID (to make sure you don't have old data)
assert st.session_state.quirky_id == current_run.quirky_id

st.write(f"######  :gray[ Analysis ID: {st.session_state.quirky_id}]")
st.title("Genomic Coordinates Mode")
st.write("Welcome to Genomic Coordinates Analysis Page. Before proceeding, "
         "please take another look at your ranges below before starting.")

# Prints on screen the ranges used for the extraction, if they are present
if current_run.regions_extracted == True:
    current_run.formatted_genomic_ranges

if current_run.regions_extracted == False:
    # Formats the ranges to use for the extraction
    try:
        current_run.formatted_genomic_ranges = format_input_df(current_run.genomic_ranges)
        current_run.formatted_genomic_ranges

    # Not the best error-handling, but will do for now...
    except (IndexError, TypeError, AttributeError):
        col1, col2 = st.columns([3, 1])
        with col1:
            st.markdown("<br><br>", unsafe_allow_html=True)
            st.warning("Please provide valid Genomic Ranges before moving to this page!")
        with col2:
            st.image("martas/error_marta.png", width=130)
        st.stop()

    if st.button("Extract Regions"):

        pre_feature_x, post_feature_y, noise_n = extract_regions_from_bw_parallel(st.session_state.bigwig_files,
                                                                                  current_run.formatted_genomic_ranges)
        # Save the three dataframes as a list, so to be able to return to the raw plots
        current_run.genomic_slices = [pre_feature_x, post_feature_y, noise_n]
        # Create a single df and coordinates, to be used by run_analysis()
        current_run.genomic_merged, current_run.genomic_merged_ranges = genomic_mode_formatter(pre_feature_x,
                                                                                                  post_feature_y,
                                                                                                  noise_n)
        current_run.regions_extracted = True

if current_run.regions_extracted == True:
    # Show Raw data
    with st.expander("Show raw count data per region"):
        st.write("""In this tab, you can explore the raw data that has been extracted for your genomic intervals of choice.
        You can expand each table by clicking on the top-right "fullscreen" button. 
        If you wish to download one of the intervals, you can also download the ranges locally as a .csv.""")
        col1, col2,col3 = st.columns(3)
        with col1:
            st.write("First Region")
            current_run.genomic_slices[0]
        with col2:
            st.write("Second Region")
            current_run.genomic_slices[1]
        with col3:
            st.write("Noise control Region")
            current_run.genomic_slices[2]

    # Show aggregated plots per genomic range extracted
    with st.expander("Show aggregate plots per region"):
        st.write("""In this tab, you can visualise the global aggregate (across all studies) for your genomic intervals of choice.
        You can expand each table by clicking on the top-right "fullscreen" button. 
        If you wish to download one of the plots, click on the '...' button on the same side.""")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.write("First Region")
            agg_df_x = global_aggregate(current_run.genomic_slices[0], True)
            # Plot
            st.line_chart(agg_df_x, x = "Position", y = "Global Aggregate")
        with col2:
            st.write("Second Region")
            agg_df_y = global_aggregate(current_run.genomic_slices[1], True)
            # Plot
            st.line_chart(agg_df_y, x = "Position", y = "Global Aggregate")
        with col3:
            st.write("Noise control Region")
            agg_df_n = global_aggregate(current_run.genomic_slices[2], True)
            # Plot
            st.line_chart(agg_df_n, x = "Position", y = "Global Aggregate")

    if st.button("Next"):
        st.switch_page("pages/page_4.py")
