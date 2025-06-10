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
else:
    st.write(f"""You selected the following ranges for the analysis:

    - First Region - From {coordinates["x_slice"][0]} to {coordinates["x_slice"][1]}
    - Second Region - From {coordinates["y_slice"][0]} to {coordinates["y_slice"][1]}
    - Noise Region - From {coordinates["n_slice"][0]} to {coordinates["n_slice"][1]}
    """)

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
    st.session_state.txid_results = run_analysis(df, coordinates,
                                            ci = st.session_state.ci,
                                            N = st.session_state.N,
                                            n_jobs = st.session_state.n_jobs)

    st.session_state.results_meta = integrate_metadata(st.session_state.txid_results,
                                                       st.session_state.metadata)

    st.session_state.analysis_ran = True


if st.session_state.get("analysis_ran", None) is True:

    with st.expander("Show Tabular Results"):

        tab1, tab2= st.tabs(["Results", "Results with metadata"])

        with tab1:
            st.session_state.txid_results
        with tab2:
            st.session_state.results_meta

    # Add expander for boxplot
    with st.expander("Interactive Plots"):

        tab3, tab4, tab5 = st.tabs(["Box plot", "Diamond plot", "Single Transcript Plot"])

        with tab3:
            fig2 = plot_rta_boxplot(st.session_state.txid_results)
            st.plotly_chart(fig2)

    # Add part for z-score plot?

    # Adds plot with breakdowns by metadata
        with tab4:

            # NB ALLOWS FOR ONLY THE FIRST 5 METADATA COLUMNS (TEMPORARY)
            meta_option = st.selectbox(
                "To visualise the diamond plot, select the medatadata column to use",
                st.session_state.results_meta.columns[:5])

            d_plot = diamond_plotter(st.session_state.results_meta, meta=meta_option)
            st.plotly_chart(d_plot, use_container_width=True)

        with tab5:
            st.write("Here you can visualise an individual Run of interest")

            # SHOULD I ADD AN OPTION TO SMOOTH THE DATA? maybe later on...

            run_option = st.multiselect("Select one or more (max 5) different Runs to plot:",
                                        st.session_state.results_meta.index)
            ml_plot = multi_line_plot(st.session_state.all_studies_df,
                                      st.session_state.utr5_length,
                                      st.session_state.cds_length,
                                      run_option)
            st.plotly_chart(ml_plot, use_container_width=True)




    #with st.expander("Diamond Plot"):
        # Insert here function to merge results with metadata...

        # Function to select by which metadata to plot the results...

        # st.button("Select"):
        # generate the plot

else:
    st.write("Here you will see the results once the analysis is complete")

# ADD FUNCTION TO EXPORT USED RANGES?