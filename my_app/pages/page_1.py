# page_1.py
from module0 import *
import streamlit as st
import pandas as pd

st.title("Genomic Coordinates Mode")

st.write("Welcome to Genomic Coordinates Analysis Page. Before proceeding, "
         "please take another look at your ranges below before starting.")

try:
    ranges = format_input_df(st.session_state.my_ranges)
    ranges

# Not the best error-handling, but will do for now
except (IndexError, TypeError, AttributeError):
    col1, col2 = st.columns([3, 1])
    with col1:
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.warning("Please provide valid Genomic Ranges before moving to this page!")
    with col2:
        st.image("martas/error_marta.png", width=130)
    st.stop()


if st.button("Extract Regions"):

    pre_feature_x, post_feature_y, noise_n = extract_regions_from_bw_parallel(st.session_state.bigwig_files, ranges)
    # Save the three dataframes as a list, so to be able to return to the raw plots
    st.session_state.genomic_slices = [pre_feature_x, post_feature_y, noise_n]
    # Create a single df and coordinates, to be used by run_analysis()
    st.session_state.genomic_merged, st.session_state.genomic_ranges = genomic_mode_formatter(pre_feature_x,
                                                                                              post_feature_y,
                                                                                              noise_n)
    st.session_state.regions_extracted = True


if "regions_extracted" in st.session_state and st.session_state.regions_extracted == True:

    # Show Raw data
    with st.expander("Show raw count data per region"):
        st.write("""In this tab, you can explore the raw data that has been extracted for your genomic intervals of choice.
        You can expand each table by clicking on the top-right "fullscreen" button. 
        If you wish to download one of the intervals, you can also download the ranges locally as a .csv.""")
        col1, col2,col3 = st.columns(3)
        with col1:
            st.write("First Region")
            st.session_state.genomic_slices[0]
        with col2:
            st.write("Second Region")
            st.session_state.genomic_slices[1]
        with col3:
            st.write("Noise control Region")
            st.session_state.genomic_slices[2]

    # Show aggregated plots per genomic range extracted
    with st.expander("Show aggregate plots per region"):
        st.write("""In this tab, you can visualise the global aggregate (across all studies) for your genomic intervals of choice.
        You can expand each table by clicking on the top-right "fullscreen" button. 
        If you wish to download one of the plots, click on the '...' button on the same side.""")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.write("First Region")
            agg_df_x = global_aggregate(st.session_state.genomic_slices[0], True)
            # Plot
            st.line_chart(agg_df_x, x = "Position", y = "Global Aggregate")
        with col2:
            st.write("Second Region")
            agg_df_y = global_aggregate(st.session_state.genomic_slices[1], True)
            # Plot
            st.line_chart(agg_df_y, x = "Position", y = "Global Aggregate")
        with col3:
            st.write("Noise control Region")
            agg_df_n = global_aggregate(st.session_state.genomic_slices[2], True)
            # Plot
            st.line_chart(agg_df_n, x = "Position", y = "Global Aggregate")

    if st.button("Next"):
        st.switch_page("pages/page_4.py")
