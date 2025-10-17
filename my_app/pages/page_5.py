# page_5.py
import streamlit as st
from module2 import *
import pandas as pd
from TranscriptID_results_helpers import *
from app_helpers import *

current_run = st.session_state.runs[st.session_state.quirky_id]
st.write(f"######  :gray[ Analysis ID: {st.session_state.quirky_id}]")
st.title("Transcript Activity Analysis")
st.write("In this page, you can run the analysis for a transcript of your choice, over the selected regions."
         "\n Please, refer to the 'Help me' expander for information on the options available.")

# Check if dataframe has been loaded
df = current_run.all_studies_df
if df is None:
    col1, col2 = st.columns([1,  3])
    with col1:
        st.image("martas/error_marta.png", width=127)
    with col2:
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.warning("Please, load the transcript before moving to this page!")
    st.stop()

# There should be a part where it checks the existence of the coordinates df
coordinates = current_run.txid_ranges
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

# WRAP ME INTO A SINGLE FUNCTION!
with st.expander("Help me"):
    st.write(" ***CI***: The confidence interval. Default is 95%.")
    st.write(" ***Number of permutations***: The number of permutations performed by the internal statistical tests."
             "A higher number of permutations allows for detection of smaller changes, but requires more time."
             "Default is 10000, while max value allowed is 100000.")
    st.write(" ***Workers to allocate***: The number of workers used by the permutation test. The greater the number, "
             "the faster the analysis. Roughly translated to cores available. Not sure if it's internally used, or"
             "if I just placed it there with the intention of implementing it later on.")


with st.expander("Change options"):
    col1, col2, col3 = st.columns([2, 2, 2])
    with col1:
        current_run.ci = st.number_input("Choose CI", min_value=1, max_value=99, value=95)
    with col2:
        current_run.N = st.number_input("Choose number of permutations", min_value=10, max_value=100000, value=10000)
    with col3:
        current_run.n_jobs = st.number_input("Choose workers to allocate", min_value=1, max_value=4, value=3)


# Sets up coordinates dictionary to be used by the function
if st.button("Run Analysis"):
    df = df.drop(columns=["Position"])

    current_run.txid_results = run_permutation_analysis(df, coordinates,
                                                        ci = current_run.ci,
                                                        N = current_run.N,
                                                        n_jobs = current_run.n_jobs)

    current_run.txid_results_meta = integrate_metadata(current_run.txid_results,
                                                       st.session_state.metadata)

    current_run.txid_analysis_done = True


if current_run.txid_analysis_done is True:

    srr_lookup_sidebar2(current_run.txid_results, current_run)

    with st.expander("Show Tabular Results"):

        tab1, tab2= st.tabs(["Detailed Results", "Results with Metadata"])

        with tab1:
            st.write("Here you can visualise the results of the analysis.")
            current_run.txid_results
        with tab2:
            st.write("Here you can visualise the results with relative metadata, either in a shortened or full version.")
            col1, col2, col3 = st.columns([2.5, 2, 6])
            verbose = False
            with col1:
                if st.button("Selected Metadata"):
                    verbose = False
            with col2:
                if st.button("Full Metadata"):
                    verbose = True

            if not verbose:    # Display fewer metadata
                st.dataframe(current_run.txid_results_meta,
                         column_order=("RTA", "padj_bh", "significant", "TISSUE", "CELL_LINE", "Info",
                                       "BioSample"))
            else:
                current_run.txid_results_meta


    # Add expander for boxplot
    with st.expander("Interactive Plots"):

        tab3, tab3_log, tab4, tab4_log,tab5, tab6 = st.tabs(
            ["Box plot (RTAs)", "Box plot (logs)",
             "Diamond plot (RTA)", "Diamond plot (logs)",
             "Single Transcript Plot", "Aggregated Plots"])

        with tab3:
            fig2 = plot_rta_boxplot(results_df = current_run.txid_results)
            st.plotly_chart(fig2)

        with tab3_log:
            # plot the results
            fig2_bis = plot_rta_boxplot(results_df = current_run.txid_results,
                                        logs =  True,
                                        box_fill = "#c44601")

            st.plotly_chart(fig2_bis,
                            key = "figure_2_bis")

    # Adds plot with breakdowns by metadata
        with tab4:

            # NB ALLOWS FOR ONLY THE FIRST 5 METADATA COLUMNS (TEMPORARY)
            meta_option = st.selectbox(
                "To visualise the diamond plot, select the medatadata column to use",
                current_run.txid_results_meta.columns[:5]) # select first 5

            d_plot = new_diamond_plotter(current_run.txid_results_meta,
                                     meta=meta_option,
                                     logs = False)
            st.plotly_chart(d_plot, use_container_width=True)

        with tab4_log:

            # NB ALLOWS FOR ONLY THE FIRST 5 METADATA COLUMNS (TEMPORARY)
            meta_option_bis = st.selectbox(
                "To visualise the diamond plot, select the medatadata column to use",
                current_run.txid_results_meta.columns[:5],
                key = "meta_dplot_bis") # select first 5

            d_plot_bis = diamond_plotter(current_run.txid_results_meta,
                                     meta=meta_option_bis,
                                     logs = True)
            st.plotly_chart(d_plot_bis,
                            use_container_width=True,
                            key = "figure_d_bis")

        with tab5:
            st.write("Here you can visualise an individual Run of interest")

            # SHOULD I ADD AN OPTION TO SMOOTH THE DATA? maybe later on...

            run_option = st.multiselect("Select one or more (max 5) different Runs to plot:",
                                        current_run.txid_results_meta.index)
            ml_plot = multi_line_plot(current_run.all_studies_df,
                                      current_run.utr5_length,
                                      current_run.cds_length,
                                      run_option,
                                      current_run.txid_ranges)
            st.plotly_chart(ml_plot, use_container_width=True)

        with tab6:
            st.write("Here you can visualise all the data aggregated (regardless of the significance of the RTA), depending on the metadata available.")

            # SHOULD I ADD AN OPTION TO SMOOTH THE DATA? maybe later on...
            st.session_state.agg_meta_option = st.selectbox(
                "Select the medatadata field to use for the aggregation",
                current_run.txid_results_meta.columns[:5])

            # Now to select the individual condition or conditions
            aggregations = st.multiselect("You can visualise one or more (max 3) different conditions at the same time",
                                        current_run.txid_results_meta[st.session_state.agg_meta_option].unique())

            # Button to choose whether to aggregate only significant ones or everything?
            agg_selection = aggregator(current_run.all_studies_df,
                                       current_run.txid_results_meta,
                                       st.session_state.agg_meta_option)

            # Function to normalise data, then aggregate it based on the elements selected in "aggregations".
            agg_meta_plot = agg_line_plot(agg_selection,
                                          current_run.utr5_length,
                                          current_run.cds_length,
                                          aggregations = aggregations,
                                          coordinates = current_run.txid_ranges)

            st.plotly_chart(agg_meta_plot, use_container_width=True)

            # Should I use the new method (smooth criminal) here?


    #with st.expander("Diamond Plot"):
        # Insert here function to merge results with metadata...

        # Function to select by which metadata to plot the results...

        # st.button("Select"):
        # generate the plot

else:
    st.write("Here you will see the results once the analysis is complete")
