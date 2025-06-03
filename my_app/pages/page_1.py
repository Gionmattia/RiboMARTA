# page_1.py
from module0 import *
import streamlit as st
import pandas as pd

st.title("Genomic Coordinates Mode")

st.write("Welcome to Genomic Coordinates Analysis Page. Before proceeding, "
         "please take another look at your ranges below before starting.")

# Load resources (should use a params.config file)
# Required Resources (should be placed into a config.file)
#bigwigs_dir = "/home/gionmattia/Desktop/OAZ3/bed_and_bigwigs/m_musculus/a_shifted_mus_musculus/a_shifted"
#bigwig_files = files = [os.path.abspath(os.path.join(bigwigs_dir, f)) for f in os.listdir(bigwigs_dir) if os.path.isfile(os.path.join(bigwigs_dir, f))]

try:
    # Add new function to format regions_df here.
    ranges = format_input_df(st.session_state.my_ranges)
    ranges

# Not the best error-handling, but will do for now
except (IndexError, TypeError, AttributeError):
    col1, col2 = st.columns([3,  1])
    with col1:
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.warning("Please provide valid Genomic Ranges before moving to this page!")
    with col2:
        st.image("martas/error_marta.png", width=130)

if st.button("Extract Regions"):

    pre_feature_x, post_feature_y, noise_n = extract_regions_from_bw_parallel(st.session_state.bigwig_files, ranges)
    st.session_state.genomic_slices = {"pre_feature": pre_feature_x, "post_feature": post_feature_y, "noise": noise_n}

    with st.expander("Show raw count data per region"):
        st.write("""In this tab, you can explore the raw data that has been extracted for your genomic intervals of choice.
        You can expand each table by clicking on the top-right "fullscreen" button. 
        If you wish to download one of the intervals, you can also download the ranges locally as a .csv.""")
        col1, col2,col3 = st.columns(3)
        with col1:
            st.write("First Region")
            pre_feature_x

        with col2:
            st.write("Second Region")
            post_feature_y

        with col3:
            st.write("Noise control Region")
            noise_n

    with st.expander("Show aggregate plots per region"):
        st.write("""In this tab, you can visualise the global aggregate (across all studies) for your genomic intervals of choice.
        You can expand each table by clicking on the top-right "fullscreen" button. 
        If you wish to download one of the plots, click on the '...' button on the same side.""")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.write("First Region")
            # aggregate
            agg_df_x = global_aggregate(pre_feature_x, True)
            # Display
            st.line_chart(agg_df_x, x = "Position", y = "Global Aggregate")

        with col2:
            st.write("Second Region")
            agg_df_y = global_aggregate(post_feature_y, True)
            # Display
            st.line_chart(agg_df_y, x = "Position", y = "Global Aggregate")

        with col3:
            st.write("Noise control Region")
            agg_df_n = global_aggregate(noise_n, True)
            # Display
            st.line_chart(agg_df_n, x = "Position", y = "Global Aggregate")