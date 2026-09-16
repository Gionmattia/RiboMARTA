# Genomic_Mode_Setup (old page_1.py)
from module0 import *
import streamlit as st
from app_helpers import *
from home_helpers import genomic_mode_formatter
from smoothing_helpers import *
from STP_helpers import *
from streamlit_bokeh import streamlit_bokeh
from Style_options import genco_mode_style

quirky_id_checkpoint()

# We are retrieving the AnalysysStatus class corresponding to the quirky_id, from the dictionary "runs"
analysis = st.session_state.analysis_runs[st.session_state.quirky_id]
# This bit tests your current id matches the AnalysisStatus ID (to make sure you don't have old data)
assert st.session_state.quirky_id == analysis.quirky_id
current_mode = analysis.genomic_mode

genco_mode_style()

st.write(f"######  :gray[ Analysis ID: {analysis.quirky_id}]")
st.title("Genomic Coordinates Mode")
st.write("Welcome to Genomic Coordinates Analysis Page. Before proceeding, "
         "please take another look at your transcript_ranges below before starting.")

# Prints on screen the transcript_ranges used for the extraction, if present
if current_mode.regions_extracted == True:
    coord_display_df = current_mode.formatted_genomic_ranges.copy()
    coord_display_df = coord_display_df.set_index("Region")
    st.table(coord_display_df)

if current_mode.regions_extracted == False:
    # Formats the transcript_ranges to use for the extraction
    try:
        current_mode.formatted_genomic_ranges = format_input_df(current_mode.genomic_ranges)
        current_mode.formatted_genomic_ranges

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

        pre_feature_x, post_feature_y, noise_n = extract_regions_from_bw_parallel(analysis.filters.selected_studies,
                                                                                  current_mode.formatted_genomic_ranges)
        # Save the three dataframes as a list, so to be able to return to the raw plots
        current_mode.genomic_slices = [pre_feature_x, post_feature_y, noise_n]

        if analysis.mappability_path is not None:
            # Transform path to mappability .bw to a list
            mapp_bw_as_list = [analysis.mappability_path]
            # Change transcript_ranges so to work with UCSC format (used by mappability track)
            ranges_formatted = current_mode.formatted_genomic_ranges.copy()
            ranges_formatted["Chromosome"] = "chr" + current_mode.formatted_genomic_ranges["Chromosome"].astype(str)
            # Extract three regions
            mapp_x, mapp_y, mapp_n = extract_regions_from_bw_parallel(mapp_bw_as_list, ranges_formatted)

            # rename column for consistency with plotting functions
            other_col = [c for c in mapp_x.columns if c not in ['Position']][0]
            mapp_x = mapp_x.rename(columns={other_col: 'Mappability'})
            mapp_y = mapp_y.rename(columns={other_col: 'Mappability'})
            mapp_n = mapp_n.rename(columns={other_col: 'Mappability'})

            # store
            current_mode.mappabilities = [mapp_x, mapp_y, mapp_n]

            # Create a single df and coordinates but for mappability, to be used by genco_data_smoother
            current_mode.genomic_merged_mappability_df, current_mode.genomic_merged_mappability_ranges = genomic_mode_formatter(
                mapp_x,
                mapp_y,
                mapp_n)

        # Create a single df and coordinates, to be used by run_analysis()
        current_mode.genomic_merged, current_mode.genomic_merged_ranges = genomic_mode_formatter(pre_feature_x,
                                                                                                  post_feature_y,
                                                                                                  noise_n)
        # Adds backup (to be used to reverse smoothing options down)
        current_mode.genomic_merged_backup = current_mode.genomic_merged.copy()

        current_mode.regions_extracted = True

if current_mode.regions_extracted == True:
    # Show Raw data
    with st.expander("Show raw count data per region"):
        st.write("""In this tab, you can explore the raw data that has been extracted for your genomic intervals of choice.
        You can expand each table by clicking on the top-right "fullscreen" button. 
        If you wish to download one of the intervals, you can also download the transcript_ranges locally as a .csv.""")
        col1, col2,col3 = st.columns(3)
        with col1:
            st.write("Baseline Region")
            current_mode.genomic_slices[0]
        with col2:
            st.write("Target Region")
            current_mode.genomic_slices[1]
        with col3:
            st.write("Noise control Region")
            current_mode.genomic_slices[2]

    # Show aggregated plots per genomic range extracted
    st.write("### Aggregate plots per region")
    st.write("""In these tabs, you can visualise the global aggregate (across all studies) for your genomic intervals of choice.
    You can expand each table by clicking on the top-right "fullscreen" button. 
    If you wish to download one of the plots, click on the '...' button on the same side.""")
    tab1, tab2, tab3 = st.tabs(["Baseline Region", "Target Region", "Noise control Region"])
    with tab1:
        st.write("Baseline Region")
        agg_df_x = global_aggregate(current_mode.genomic_slices[0],
                                    normalize = False,
                                    keep_pos = True)

        if current_mode.mappabilities is not None:

            x_plot = colorless_transcript_plot(df = agg_df_x,
                                               y_axis = "Global Aggregate",
                                               mappability_df = current_mode.mappabilities[0])
            streamlit_bokeh(x_plot, use_container_width=False)

        else:
            x_plot = colorless_transcript_plot(df = agg_df_x,
                                               y_axis = "Global Aggregate")
            streamlit_bokeh(x_plot, use_container_width=False)

    with tab2:
        st.write("Target Region")
        agg_df_y = global_aggregate(current_mode.genomic_slices[1],
                                    normalize = False,
                                    keep_pos = True)

        if current_mode.mappabilities is not None:

            y_plot = colorless_transcript_plot(df = agg_df_y,
                                               y_axis = "Global Aggregate",
                                               mappability_df = current_mode.mappabilities[1])
            streamlit_bokeh(y_plot, use_container_width=False)

        else:
            y_plot = colorless_transcript_plot(df = agg_df_y,
                                               y_axis = "Global Aggregate")
            streamlit_bokeh(y_plot, use_container_width=False)

    with tab3:
        st.write("Noise control Region")
        agg_df_n = global_aggregate(current_mode.genomic_slices[2],
                                    normalize = False,
                                    keep_pos = True)

        if current_mode.mappabilities is not None:

            n_plot = colorless_transcript_plot(df = agg_df_n,
                                               y_axis = "Global Aggregate",
                                               mappability_df = current_mode.mappabilities[2])
            streamlit_bokeh(n_plot, use_container_width=False)

        else:
            n_plot = colorless_transcript_plot(df = agg_df_n,
                                               y_axis = "Global Aggregate")
            streamlit_bokeh(n_plot, use_container_width=False)


    with st.expander("Select Smoothing Option"):

        # Failsafe - If no mappability is given, no smoothing is done
        if current_mode.mappabilities is None:
            valid_options = ["No Smoothing", "Criminal Smoothing (Experimental)"]
        else:
            valid_options = ["No Smoothing",
                     "Mappability Smoothing (Proportional)",
                     "Mappability Smoothing (Aggressive)",
                     "Criminal Smoothing (Experimental)"]
        choice = st.radio(
            "_",
            label_visibility="collapsed",
            options= valid_options,
            horizontal=True,  # makes it look more toggle-like
            captions=["Default option. No smoothing is applied.",
                      "Smoothes the coverage per position by their mappability value (see plot).",
                      "Sets all non-1 mappability positions to 0-counts.",
                      "Per transcript region, finds extremes and smoothes them. Z-score 4 is a good start."],
        )


        if choice == "Criminal Smoothing (Experimental)":
            z_score_threshold = st.slider("z_score_threshold", 1, 10)
            c1, c2 = st.columns(2)
            with c1:
                st.caption("Abrasive (removes a lot)")
            with c2:
                st.caption("<div style='text-align:right'>Permissive (allows peaks) </div>", unsafe_allow_html=True)

        if st.button("Smooth it!"):
            current_mode.smoothing_option = choice
            if choice == "Criminal Smoothing (Experimental)":
                current_mode.smoothing_criminal_thresh = z_score_threshold
            # overwrites all_studies_df (ndr backup still safe) to add smoothing options
            genco_data_smoother(current_mode, current_mode.smoothing_option)


            smoothed_w_ga = global_aggregate(current_mode.genomic_merged,
                                        normalize=False,
                                        keep_pos=True)

            x_boundaries = current_mode.genomic_merged_ranges["x_slice"]
            x_smooth = smoothed_w_ga[x_boundaries[0]:x_boundaries[1]]

            tab4, tab5, tab6 = st.tabs(["Baseline Region", "Target Region", "Noise control Region"])

            with tab4:

                x_boundaries = current_mode.genomic_merged_ranges["x_slice"]
                x_smooth = smoothed_w_ga[x_boundaries[0]:x_boundaries[1]]

                after_smooth_plot_x = colorless_transcript_plot(df = x_smooth,
                                                                y_axis = "Global Aggregate",
                                                                mappability_df = current_mode.mappabilities[0])

                streamlit_bokeh(after_smooth_plot_x, use_container_width=False)

            with tab5:

                y_boundaries = current_mode.genomic_merged_ranges["y_slice"]
                y_smooth = smoothed_w_ga[y_boundaries[0]:y_boundaries[1]]

                after_smooth_plot_y = colorless_transcript_plot(df = y_smooth,
                                                                y_axis = "Global Aggregate",
                                                                mappability_df = current_mode.mappabilities[1])

                streamlit_bokeh(after_smooth_plot_y, use_container_width=False)


            with tab6:

                n_boundaries = current_mode.genomic_merged_ranges["n_slice"]
                n_smooth = smoothed_w_ga[n_boundaries[0]:n_boundaries[1]]

                n_smooth

                after_smooth_plot_n = colorless_transcript_plot(df = n_smooth,
                                                                y_axis = "Global Aggregate",
                                                                mappability_df = current_mode.mappabilities[2])

                streamlit_bokeh(after_smooth_plot_n, use_container_width=False)

            # Prepare data
            # ga_df = global_aggregate(all_studies_df, normalize=normalize)
            # processed_df, k = plot_smoother(ga_df, smoothing)
            # #


    if st.button("Next"):
        st.switch_page("pages/Genomic_Mode_Results.py")
