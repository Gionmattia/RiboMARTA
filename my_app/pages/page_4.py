# page_4.py
import streamlit as st
from module2 import *
import pandas as pd


st.title("Transcript Activity Analysis")

st.write("In this page, you can run the analysis for a transcript of your choice, over the selected regions.")

# Check if dataframe has been loaded
df = st.session_state.get("all_studies_df", None)
if df is None:
    col1, col2 = st.columns([1,  3])
    with col1:
        st.image("martas/error_marta.png", width=127)
    with col2:
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.warning("Please, load the transcript before moving to this page!")
    st.stop()


# There should be a part where it checks the existence of the coordinates df
coordinates = st.session_state.get("coordinates", None)
if coordinates is None:
    col1, col2 = st.columns([3,  1])
    with col1:
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.warning("Please, select ranges of transcript to compare before proceeding!")
    with col2:
        st.image("martas/error_marta.png", width=127)
    st.stop()

with st.expander("Change options"):
    col1, col2, col3 = st.columns([2, 2, 2])
    with col1:
        st.session_state.ci = st.number_input("Choose CI", min_value=1, max_value=99, value=95)
    with col2:
        st.session_state.N = st.number_input("Choose number permutations", min_value=10, max_value=100000, value=10000)
    with col3:
        st.session_state.n_jobs = st.number_input("Choose workers to allocate", min_value=1, max_value=4, value=3)

# Sets up coordinates dictionary to be used by the function
if st.button("Run Analysis"):
    df = df.drop(columns=["Position"])
    st.session_state.results = run_analysis(df, coordinates,
                                            ci = st.session_state.ci,
                                            N = st.session_state.N,
                                            n_jobs = st.session_state.n_jobs)

    with st.expander("Show results"):
        st.session_state.results

    # Add expander for boxplot
    with st.expander("Interactive boxplot"):
        fig2 = plot_observed_rta_interactive(st.session_state.results)
        st.plotly_chart(fig2)
    # Add part for z-score plot?

# ADD FUNCTION TO EXPORT USED RANGES?