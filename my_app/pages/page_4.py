# page_4.py
import streamlit as st
from module2 import *
import pandas as pd
from TranscriptID_results_helpers import *
from genomic_mode_helpers import *

# We are retrieving the AnalysysStatus class corresponding to the quirky_id, from the dictionary "runs"
analysis = st.session_state.analysis_runs[st.session_state.quirky_id]
# This bit tests your current id matches the AnalysisStatus ID (to make sure you don't have old data)
current_mode = analysis.genomic_mode
st.write(f"######  :gray[ Analysis ID: {analysis.quirky_id}]")

st.title("Genomic Coodrinates Activity Analysis")
st.write("In this page, you can run the analysis for your selected genomic coordinates of choice.")

# Check if dataframe has been loaded
df = current_mode.genomic_merged
genomic_upload = analysis.genomic_upload
if df is None:
    col1, col2 = st.columns([1,  3])
    with col1:
        st.image("martas/error_marta.png", width=127)
    with col2:
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.warning("Please, select the genomic ranges and extract them before moving to this page!")
    st.stop()


# There should be a part where it checks the existence of the coordinates df
coordinates =  current_mode.genomic_merged_ranges
if coordinates is None:
    col1, col2 = st.columns([3,  1])
    with col1:
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.warning("Error. Genomic ranges still not generated!")
    with col2:
        st.image("martas/error_marta.png", width=127)
elif genomic_upload == True:
    st.write("You are currently in: Genomic Coordinates - File Upload mode")

else:
    st.write(f"""You selected the following genomic coordinates for the analysis:

    - First Region - { current_mode.genomic_ranges[0] }
    - Second Region - { current_mode.genomic_ranges[1] }
    - Noise Region - { current_mode.genomic_ranges[2] }
    """)

with st.expander("Change options"):
    col1, col1_bis, col2, col3, col4, col5 = st.columns([2, 2, 2, 2, 2, 2])
    with col1:
        current_mode.alpha_baseline = st.number_input("Choose alpha for baseline", min_value=0.00, max_value=0.99, value=0.05)
    with col1_bis:
        current_mode.alpha_target = st.number_input("Choose alpha for RTA test", min_value=0.0, max_value=0.99, value=0.05)
    with col2:
        current_mode.ci = st.number_input("Choose CI", min_value=1, max_value=99, value=95)
    with col3:
        current_mode.N = st.number_input("Choose # of permutations", min_value=10, max_value=100000, value=10000)
    with col4:
        current_mode.n_jobs = st.number_input("Choose workers to allocate", min_value=1, max_value=4, value=3)
    with col5:
        use_seed = st.checkbox("Set random seed for reproducibility", value=False)
        if use_seed:
            current_mode.seed = st.number_input("Random seed", min_value=0, value=42)
        else:
            current_mode.seed = None

# Sets up coordinates dictionary to be used by the function
if st.button("Run Analysis"):
    # NEED TO SAVE TO THE ANALYSIS STATUS THE ANALYSIS VALUES


    df = df.drop(columns=["Position"])
    current_mode.genco_results, current_mode.seed = run_permutation_analysis(df = df, coordinates = coordinates,
                                                                            alpha_baseline = current_mode.alpha_baseline,
                                                                            alpha_target = current_mode.alpha_target,
                                                                            ci = current_mode.ci,
                                                                            N = current_mode.N,
                                                                            n_jobs = current_mode.n_jobs,
                                                                            random_seed = current_mode.seed)

    current_mode.results_meta = integrate_metadata(current_mode.genco_results, # DOUBLE CHECK THIS
                                                    st.session_state.metadata)

    current_mode.genco_analysis_done = True


if current_mode.genco_analysis_done == True:

    with st.expander("Show Tabular Results"):

        tab1, tab2= st.tabs(["Detailed Results", "Results with Metadata"])

        with tab1:
            st.write("Here you can visualise the results of the analysis.")
            current_mode.genco_results
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
                st.dataframe(current_mode.results_meta,
                         column_order=("RTA", "padj_bh", "significant", "TISSUE", "CELL_LINE", "Info",
                                       "BioSample"))
            else:
                current_mode.results_meta


    # Add expander for boxplot
    with st.expander("Interactive Plots"):

        tab3, tab3_log, tab4, tab4_log, tab5, tab6 = st.tabs(
            ["Box plot (RTAs)", "Box plot (logs)",
             "Diamond plot (RTA)", "Diamond plot (logs)",
             "Single Transcript Plot", "Aggregated Plots"])

        with tab3:
            fig2 = plot_rta_boxplot(current_mode.genco_results)
            st.plotly_chart(fig2)

            if st.button("Store in the report", key = "rta_boxplot"):
                analysis.user_plots.append(fig2)
                st.success(f"Plot added! You currently have {len(analysis.user_plots)} plots.")


        with tab3_log:
            # Modify results so that they are log2-based now!
            log_2_genco_results = current_mode.genco_results
            # plot the results
            fig2_bis = plot_rta_boxplot(log_2_genco_results, logs = True, box_fill = "#c44601")
            st.plotly_chart(fig2_bis, key = "plot_bis")

            if st.button("Store in the report", key = "log_rta_boxplot"):
                analysis.user_plots.append(fig2_bis)
                st.success(f"Plot added! You currently have {len(analysis.user_plots)} plots.")


    # Adds plot with breakdowns by metadata
        with tab4:
            # NB ALLOWS FOR ONLY THE FIRST 5 METADATA COLUMNS (TEMPORARY)
            meta_option = st.selectbox(
                "To visualise the diamond plot, select the medatadata column to use",
                current_mode.results_meta.columns[:5]) # select first 5

            d_plot = new_diamond_plotter(current_mode.results_meta,
                                         meta=meta_option,
                                         logs = False)
            st.plotly_chart(d_plot, use_container_width=True)

            if st.button("Store in the report", key = "diamond plot"):
                analysis.user_plots.append(d_plot)
                st.success(f"Plot added! You currently have {len(analysis.user_plots)} plots.")

        with tab4_log:

            # NB ALLOWS FOR ONLY THE FIRST 5 METADATA COLUMNS (TEMPORARY)
            meta_option_bis = st.selectbox(
                "To visualise the diamond plot, select the medatadata column to use",
                current_mode.results_meta.columns[:5],
                key = "meta_dplot_bis") # select first 5

            d_plot_bis = diamond_plotter(current_mode.results_meta,
                                     meta=meta_option_bis,
                                     logs = True)
            st.plotly_chart(d_plot_bis,
                            use_container_width=True,
                            key = "figure_d_bis")

            if st.button("Store in the report", key = "log diamond plot"):
                analysis.user_plots.append(d_plot_bis)
                st.success(f"Plot added! You currently have {len(analysis.user_plots)} plots.")

        with tab5:
            st.write("Here you can visualise an individual Run of interest")

            run_option = st.multiselect("Select one or more (max 5) different Runs to plot:",
                                        current_mode.results_meta.index)

            ml_plot = multi_line_plot_genomic(current_mode.genomic_merged,
                                      run_option,
                                      coordinates)
            st.plotly_chart(ml_plot, use_container_width=True)

        with tab6:
            st.write("Here you can visualise all the data aggregated (regardless of the significance of the rta), depending on the metadata available.")

            # SHOULD I ADD AN OPTION TO SMOOTH THE DATA? maybe later on...
            agg_meta_option = st.selectbox(
                "Select the medatadata field to use for the aggregation",
                current_mode.results_meta.columns[:5])

            # SAME AS ABOVE< THIS SHOULD BE MADE TO WORK WITH THREE DFS, NOT JUST ONE!
            # # Now to select the individual condition or conditions
            aggregations = st.multiselect("You can visualise one or more (max 3) different conditions at the same time",
                                         current_mode.results_meta[agg_meta_option].unique())

            # Button to choose whether to aggregate only significant ones or everything.
            agg_selection = aggregator(current_mode.genomic_merged,
                                        current_mode.results_meta,
                                        agg_meta_option)

             # Function to normalise data, then aggregate it based on the elements selected in "aggregations".
            genom_agg_meta_plot = genom_agg_line_plot(agg_selection,
                                          aggregations = aggregations,
                                          coordinates = coordinates)

            st.plotly_chart(genom_agg_meta_plot, use_container_width=True)

            # Should I use the new method (smooth criminal) here?


else:
    st.write("Here you will see the results once the analysis is complete")
