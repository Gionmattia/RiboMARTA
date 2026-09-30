# Genomic_mode_results.py (old page_4.py)
import streamlit as st
from module2 import *
from TranscriptID_results_helpers import *
from genomic_mode_helpers import *
from Style_options import genco_mode_style
from streamlit_bokeh import streamlit_bokeh

from STP_helpers import *
/ from bokeh.plotting import figure
from bokeh.models import (ColumnDataSource, HoverTool, PanTool, Range1d, Span, Label, Legend
)
import pandas as pd
# from bokeh.models import LinearAxis
from app_helpers import AnalysisStatus, srr_genco_lookup_sidebar, quirky_id_checkpoint, graceful_failure
# from bokeh.layouts import column
# from bokeh.models import BoxAnnotation

genco_mode_style()

quirky_id_checkpoint()

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
        st.image(st.session_state["marta_images"]["error"], width=127)
    with col2:
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.warning("Please, select the genomic transcript_ranges and extract them before moving to this page!")
    st.stop()


# There should be a part where it checks the existence of the coordinates df
coordinates =  current_mode.genomic_merged_ranges
if coordinates is None:
    col1, col2 = st.columns([3,  1])
    with col1:
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.warning("Error. Genomic transcript_ranges still not generated!")
    with col2:
        st.image(st.session_state["marta_images"]["error"], width=127)
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

    # Safety Check before execution
    bstart, bend = coordinates["x_slice"]
    baseline = df.iloc[bstart:bend]
    if (baseline.mean(axis=1) == 0).all():
        text = "There is no data from the Baseline region across all studies!"
        graceful_failure(text)

    current_mode.genco_results, current_mode.seed = run_permutation_analysis(df = df, coordinates = coordinates,
                                                                            alpha_baseline = current_mode.alpha_baseline,
                                                                            alpha_target = current_mode.alpha_target,
                                                                            ci = current_mode.ci,
                                                                            N = current_mode.N,
                                                                            n_jobs = current_mode.n_jobs,
                                                                            random_state = current_mode.seed)

    current_mode.results_meta = integrate_metadata(current_mode.genco_results, # DOUBLE CHECK THIS
                                                    st.session_state.metadata)

    current_mode.genco_analysis_done = True


if current_mode.genco_analysis_done == True:

    srr_genco_lookup_sidebar(current_mode.results_meta, current_mode, analysis)

    with st.expander("Show Tabular Results"):

        tab1, tab2, tab3= st.tabs(["Analysis Results", "Results with Selected Metadata", "Result with Full Metadtaa"])

        with tab1:
            st.write("Here you can visualise the results of the analysis. This table is useful to check "
                     "the numerical aspects of the analysis results.")
            current_mode.genco_results
        with tab2:
            st.write("Here you can visualise the results with selected metadata of your choice. Useful to have a more 'easy-to-parse' table."
                     "By default, it always showcases the SRR, the RTA, the padj, and the Total Counts.")
            all_meta_options = current_mode.results_meta.columns
            meta_options = st.multiselect("Select metadata field to visualise", all_meta_options)
            default_meta = ["RTA", "padj_bh", "Tot Counts"]
            table_meta = default_meta + meta_options
            st.dataframe(current_mode.results_meta,
                         column_order=(table_meta))

        with tab3:
            st.write("Here you can visualise the results of the analysis with the complete metadata. Quite verbose but the most informative")
            current_mode.results_meta


    # Add expander for boxplot
    with st.expander("Interactive Plots"):

        max_val = max(int(max(current_mode.genco_results["Tot Counts"])),1)


        tab3, tab4, tab5 = st.tabs(
            ["Box plot (RTAs)", "Diamond plot (RTA)",
             "Single Transcript Plot" ])

        with tab3:
            plot_logs = False
            box_colour = "#00bf7d"

            col1, col2, col3 = st.columns([2, 2, 6])
            with col1:
                if st.checkbox("Plot Log(2) values"):
                    plot_logs = True
                    box_colour = "#c44601"

                if st.checkbox("Plot significant RTAs only"):
                    data = current_mode.genco_results[current_mode.genco_results["significant"] == True]
                else:
                    data = current_mode.genco_results


            with col3:

                tot_counts_min = st.slider("Select min Tot Counts for a study to be considered",
                          0, max_val, key =  "boxplot_slider")

            data = data[data["Tot Counts"] > tot_counts_min]

            fig2 = plot_rta_boxplot(results_df = data,
                                    box_fill = box_colour,
                                    logs = plot_logs)
            st.plotly_chart(fig2)

            if st.button("Store in the report", key = "rta_boxplot"):
                analysis.user_plots.append(fig2)
                st.success(f"Plot added! You currently have {len(analysis.user_plots)} plots.")


    # Adds plot with breakdowns by metadata
        with tab4:

            if current_mode.aggregation_option == "Aggregate by metadata field":
                meta_option = "Group"
            elif current_mode.aggregation_option == "Total Aggregation":
                meta_option = "Group"
            else:
                meta_option = st.selectbox("To visualise the diamond plot, select the medatadata column to use",
                                           current_mode.results_meta.columns[:])

            show_nonsignificant =  False
            plot_logs_d = False

            col1, col2, col3 = st.columns([2, 2, 6])
            with col1:
                if st.checkbox("Plot Non Significant"):
                    show_nonsignificant = True

                if st.checkbox("Plot Log(2) values", key = "diamond_plot_"):
                    plot_logs_d = True

            with col3:

                tot_counts_min = st.slider("Select min Tot Counts for a study to be considered",
                          0, max_val, key =  "diamond_plot_slider")

            data_diamond_plot = current_mode.results_meta
            data_diamond_plot = data_diamond_plot[data_diamond_plot["Tot Counts"] > tot_counts_min]
            data_diamond_plot = data_diamond_plot.sort_values(by = meta_option, ascending = True)


            if current_mode.aggregation_option is None or current_mode.aggregation_option == "":
                jitter = 0.3
                show_n_members =  False
            else:
                jitter = 0
                show_n_members = True
                # here we add the info on the # elements per group
                data_diamond_plot = pd.merge(data_diamond_plot, current_mode.members_per_group, on = meta_option)

            p = new_diamond_plotter_bokeh(df = data_diamond_plot,
                                          meta=meta_option,
                                          logs=plot_logs_d,
                                          jitter_amount= jitter,
                                          show_nonsignificant=show_nonsignificant,
                                          show_n_members = show_n_members)
            streamlit_bokeh(p, use_container_width=True, theme="light_minimal")


            if st.button("Store in the report", key = "diamond plot"):
                analysis.user_plots.append(p)
                st.success(f"Plot added! You currently have {len(analysis.user_plots)} plots.")


        with tab5:
            st.write("Here you can visualise an individual Run of interest, in RAW COUNTS")

            tot_datasets = set(df.columns)
            if "Position" in tot_datasets:
                tot_datasets.remove("Position")
            baseline_passed = set(current_mode.results_meta.index)
            baseline_failed = tot_datasets - baseline_passed

            options_STP  = st.radio(
                "_",
                label_visibility="collapsed",
                options=["All datasets", "Baseline Expressed","Baseline Not Expressed"],
                horizontal=True,  # makes it look more toggle-like
                captions=[f"({len(tot_datasets)} datasets)",
                          f"Only datasets where Signal Baseline > Noise Control ({len(baseline_passed)})",
                          f"Only datasets where Signal Baseline indistinguishable from Noise Control ({len(baseline_failed)})",],
            )

            if options_STP == "All datasets":
                dat_selected = tot_datasets
            elif options_STP == "Baseline Expressed":
                dat_selected = baseline_passed
            elif options_STP == "Baseline Not Expressed":
                dat_selected = baseline_failed

            run_option = st.multiselect("Select one or more (max 5) different Runs to plot:",
                                        dat_selected,
                                        key = "new_multiline_stp_key")

            val = int(df.index.max()) - int(df.index.min())

            margin = max(1, val) * 0.01

            genco_ml_plot = genco_single_transcript_plot(df = df,
                                                         x_range = Range1d(start=int(df.index.min()) - margin, end=int(df.index.max()) + margin),
                                                         run_option = run_option,
                                                         mappability_df = None,
                                                         genomic_merged_ranges = current_mode.genomic_merged_ranges
                                                         )

            streamlit_bokeh(genco_ml_plot, use_container_width=False)

            if st.button("Store in the report", key = "genco_ml_plot"):
                analysis.user_plots.append(genco_ml_plot)
                st.success(f"Plot added! You currently have {len(analysis.user_plots)} plots.")


else:
    st.write("Here you will see the results once the analysis is complete")
