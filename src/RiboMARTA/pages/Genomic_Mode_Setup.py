# Genomic_Mode_Setup (old page_1.py)
from module0 import *
import streamlit as st
from app_helpers import *
from home_helpers import genomic_mode_formatter
from smoothing_helpers import *
from STP_helpers import *
from streamlit_bokeh import streamlit_bokeh
from Style_options import genco_mode_style
from genomic_mode_helpers import get_sequence

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
            st.image(st.session_state["marta_images"]["error"], width=130)
        st.stop()

    if st.button("Extract Regions"):

        baseline_region, target_region, noise_region = extract_regions_from_bw_parallel(analysis.filters.selected_studies,
                                                                                        current_mode.formatted_genomic_ranges)
        # Save the three dataframes as a list, so to be able to return to the raw plots
        current_mode.genomic_slices = [baseline_region, target_region, noise_region]

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
        current_mode.genomic_merged, current_mode.genomic_merged_ranges = genomic_mode_formatter(baseline_region,
                                                                                                 target_region,
                                                                                                 noise_region)
        # Adds backup (to be used to reverse smoothing options down)
        current_mode.genomic_merged_backup = current_mode.genomic_merged.copy()

        # Extract and store the genomic sequence for each region
        fasta = st.session_state.fasta
        sequence_dict = get_sequence(fasta=fasta, df=current_mode.formatted_genomic_ranges)
        current_mode.regions_sequence = sequence_dict

        # Insert sequences into the dataframes
        current_mode.genomic_slices[0].insert(0, 'Sequence', list(current_mode.regions_sequence["Baseline Region"]))
        current_mode.genomic_slices[1].insert(0, 'Sequence', list(current_mode.regions_sequence["Target Region"]))
        current_mode.genomic_slices[2].insert(0, 'Sequence', list(current_mode.regions_sequence["Noise Region"]))

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
            # Add download button?
        with col2:
            st.write("Target Region")
            current_mode.genomic_slices[1]
            # Add download button?
        with col3:
            st.write("Noise Control Region")
            current_mode.genomic_slices[2]
            # Add download button?

    # Show aggregated plots per genomic range extracted
    st.write("### Aggregate plots per region")
    st.write("""In these tabs, you can visualise the global aggregate (across all studies) for your genomic intervals of choice.
    You can expand each table by clicking on the top-right "fullscreen" button. 
    If you wish to download one of the plots, click on the '...' button on the same side.""")
    st.markdown("#### IMPORTANT: ")
    st.markdown("Pseudo-ORF plots were added. For each region, frames are defined from the start of each interval."
                "This is done ONLY to allow users to see the presence of in-frame starts or stops codons WITHIN a region."
                " THESE ARE USER-DEFINED GENOMIC REGIONS, THE 'READING FRAMES' ARE JUST AN ABSTRACTION TO VISUALISE BETTER!")

    tab1, tab2, tab3 = st.tabs(["Baseline Region", "Target Region", "Noise Control Region"])
    with tab1:
        st.write("Baseline Region")
        smooth = st.toggle("Smooth data", key = "baseline_toggle")
        baseline_plot = plot_ga_with_frames1(region_df = current_mode.genomic_slices[0],
                                             sequence = current_mode.regions_sequence["Baseline Region"],
                                             smoothing = smooth,
                                             mappability = current_mode.mappabilities[0].reset_index(names="X")
                                             )

        streamlit_bokeh(baseline_plot, use_container_width=False)

        if st.button("Store in the report", key="genco_baseline_plot"):
            analysis.user_plots.append(baseline_plot)
            st.success(f"Plot added! You currently have {len(analysis.user_plots)} plots.")


    with tab2:
        st.write("Target Region")
        smooth = st.toggle("Smooth data", key = "target_toggle")
        target_plot = plot_ga_with_frames1(region_df = current_mode.genomic_slices[1],
                                           sequence = current_mode.regions_sequence["Target Region"],
                                           smoothing = smooth,
                                           mappability=current_mode.mappabilities[1].reset_index(names="X")
                                           )

        streamlit_bokeh(target_plot, use_container_width=False)

        if st.button("Store in the report", key="genco_target_plot"):
            analysis.user_plots.append(target_plot)
            st.success(f"Plot added! You currently have {len(analysis.user_plots)} plots.")

    with tab3:
        st.write("Noise Control Region")
        smooth = st.toggle("Smooth data", key = "noise_toggle")
        noise_plot = plot_ga_with_frames1(region_df = current_mode.genomic_slices[2],
                                          sequence = current_mode.regions_sequence["Noise Region"],
                                          smoothing = smooth,
                                          mappability=current_mode.mappabilities[2].reset_index(names="X")
                                          )

        streamlit_bokeh(noise_plot, use_container_width=False)

        if st.button("Store in the report", key="genco_noise_plot"):
            analysis.user_plots.append(noise_plot)
            st.success(f"Plot added! You currently have {len(analysis.user_plots)} plots.")


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

            tab4, tab5, tab6 = st.tabs(["Baseline Region", "Target Region", "Noise Control Region"])

            with tab4:

                x_boundaries = current_mode.genomic_merged_ranges["x_slice"]

                df1 = current_mode.genomic_merged[x_boundaries[0]:x_boundaries[1]].copy()
                df1 = df1.drop(columns = ["Position"]).reset_index(drop = True)
                smooth_baseline_plot = plot_ga_with_frames1(region_df=df1,
                                                     sequence=current_mode.regions_sequence["Baseline Region"],
                                                     smoothing=False,
                                                     mappability=current_mode.mappabilities[0].reset_index(names="X")
                                                     )

                streamlit_bokeh(smooth_baseline_plot, use_container_width=False)


            with tab5:

                y_boundaries = current_mode.genomic_merged_ranges["y_slice"]
                df2 = current_mode.genomic_merged[y_boundaries[0]:y_boundaries[1]].copy()
                df2 = df2.drop(columns = ["Position"]).reset_index(drop = True)
                smooth_target_plot = plot_ga_with_frames1(region_df= df2,
                                                   sequence=current_mode.regions_sequence["Target Region"],
                                                   smoothing=False,
                                                   mappability=current_mode.mappabilities[1].reset_index(names="X")
                                                   )

                streamlit_bokeh(smooth_target_plot, use_container_width=False)



            with tab6:

                n_boundaries = current_mode.genomic_merged_ranges["n_slice"]
                df3 = current_mode.genomic_merged[n_boundaries[0]:n_boundaries[1]].copy()
                df3 = df3.drop(columns = ["Position"]).reset_index(drop = True)
                smooth_noise_plot = plot_ga_with_frames1(region_df=df3,
                                                  sequence=current_mode.regions_sequence["Noise Region"],
                                                  smoothing=False,
                                                  mappability=current_mode.mappabilities[2].reset_index(names="X")
                                                  )

                streamlit_bokeh(smooth_noise_plot, use_container_width=False)


    if st.button("Next"):
        st.switch_page("pages/Genomic_Mode_Results.py")
