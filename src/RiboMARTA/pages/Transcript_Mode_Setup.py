# page_2.py
import streamlit as st
from smoothing_helpers import *
from app_helpers import *
from module0 import *
from transcriptID_mode_helpers import *
from transcriptID_mode_helpers2 import *
from STP_helpers import *
from streamlit_bokeh import streamlit_bokeh
from Style_options import txid_mode_style, bubble
import streamlit as st

##### Imports CSS to style the whole page #####
txid_mode_style()

##### Analysis Checkpoints and Setup #####
#First checkpoint. Stops rendering if the user has not initialised any analysis in the Home page
quirky_id_checkpoint()

# We are retrieving the AnalysysStatus class corresponding to the quirky_id, from the dictionary "runs"
analysis = st.session_state.analysis_runs[st.session_state.quirky_id]
# This bit tests your current id matches the AnalysisStatus ID (to make sure you don't have old data)
assert st.session_state.quirky_id == analysis.quirky_id
current_mode = analysis.transcript_mode
# Retrieve the AnalysysStatus class corresponding to the quirky_id, from the dictionary "runs"
fasta = st.session_state.fasta
st.write(f"######  :gray[ Analysis ID: {st.session_state.quirky_id}]")
st.title("Transcript ID Mode")


##### Analysis Start #####

# Renders the info message for this mode
txid_mode_intro()
# Checkpoint - checks for presence of a txid (eg. ENST...) if user has not submitted one, returns error.
txid_checkpoint(current_mode)
# Checkpoint - Appears only if not in Upload mode. Button to start extraction of data.
if analysis.tx_upload == False:
    if st.button("Load Transcript"):
        current_mode.loading_tx_cov = True  # Logical flag

## All the code to extract info - needs some order ##
# Extract the data coverage info, transcript sequence, and mappability. Includes loader.
if current_mode.loading_tx_cov:

    with st.spinner("I am retrieving the coverage data..."):
        # Internal Check. Stops execution if the txid is not found for any reason.
        try:
            # 1 Retrieve transcript transcript_ranges
            transcript_ranges = extract_tx_pyranges(st.session_state.gr, current_mode.txid)

        except KeyError:
            st.error("I am sorry, but the Transcript ID given was not found in the annotation file."
                     "Please, restart from the Home page.")
            st.stop()

        # 2  Retrieve transcript features (CDS, UTRs) transcript_ranges and strandedness
        all_ranges, strandedness = regions_df_formatter(transcript_ranges)
        st.session_state.strandedness = strandedness

        # 3  Retrieve features (CDS, UTRs) lengths
        current_mode.utr5_length, current_mode.cds_length, current_mode.utr3_length = extract_features_length(
            transcript_ranges, strandedness)

        # 4 Extract coverage from given bigwigs
        try:
            current_mode.all_studies_df = extract_from_multiple_bigwigs(analysis.filters.selected_studies,
                                                                       all_ranges,
                                                                       strandedness,
                                                                       n_jobs=1)
            # Updates logic flags within current_run
            current_mode.loading_tx_cov = False
            current_mode.tx_cov_extracted = True
        except:
            st.error("I was unable to retrieve the information from the .BigWig files")
            st.stop()

    # Extracts mappability info, if available
    if analysis.mappability_path is not None:
        mappability_df = txid_map_extraction(current_mode.txid, analysis.mappability_path, st.session_state.gr)
        mappability_df["X"] = mappability_df.index
        other_col = [c for c in mappability_df.columns if c not in ['Position', 'X']][0]
        current_mode.mappability_df = mappability_df.rename(columns={other_col: 'Mappability'})

    current_mode.aggregation = None
    current_mode.smoothing = False

    # Retrieves the sequence
    current_mode.transcript_sequence = get_sequence(fasta, all_ranges, strandedness)

    # Creates back-ups for the smoothing options and aggregation options.
    current_mode.all_studies_df_backup = current_mode.all_studies_df.copy()
    current_mode.all_studies_df_pre_aggregation_backup = current_mode.all_studies_df.copy()


# Checkpoint - Checks whether the coverage has been extracted
if current_mode.tx_cov_extracted:
    # If data has been extracted, allows for download in .csv format
    if analysis.tx_upload == False:
        format_txid_df_download(current_mode)
    # If it's "File Upload Mode", reminds it to the user.
    elif analysis.tx_upload == True:
        st.markdown("### You are currently in: Ensembl Transcript ID - File Upload mode")
else:
    st.stop()

#####  Display the data once extracted #####
if current_mode.tx_cov_extracted == True and current_mode.transcript_sequence != "":
    if current_mode.triplet_ga_plot is None:
        current_mode.triplet_ga_plot = plot_ga_with_frames2(current_mode,
                                                            current_mode.all_studies_df,
                                                            y_column = "Global Aggregate smoothed",
                                                            normalize=False)

    if current_mode.triplet_ga_plot_raw is None:
        current_mode.triplet_ga_plot_raw = plot_ga_with_frames2(current_mode,
                                                            current_mode.all_studies_df,
                                                            smoothing = False,
                                                            y_column = "Global Aggregate",
                                                            normalize=False)

    if current_mode.norm_triplet_ga_plot is None:
        current_mode.norm_triplet_ga_plot = plot_ga_with_frames2(current_mode,
                                                                 current_mode.all_studies_df,
                                                                 y_column = "Global Aggregate smoothed",
                                                                 normalize=True)

    # Expanders to render the plots, if present
    ga_plots_UI_render(current_mode, analysis)


    #####  Ranges Selection #####
    with st.container(border = True):
        st.write("#### Select Ranges")
        col1, col2, col3, col4 = st.columns([4, 4, 4, 8])
        full_mRNA_length = current_mode.utr5_length + current_mode.cds_length + current_mode.utr3_length

        with col1:
            st.write("##### Baseline Region")
            st.write("\n")
            current_mode.x_slice = ranges_selector(full_mRNA_length, "x")

        with col2:
            st.write("##### Target Region")
            st.write("\n")
            current_mode.y_slice = ranges_selector(full_mRNA_length, "y")

        with col3:
            st.write("##### Noise Control Region")
            st.write("\n")
            current_mode.n_slice = ranges_selector(full_mRNA_length, "n")

        with col4:
            st.write("##### Guidelines:")
            st.write(""" When possible, pick a region 25nt away from start/stop sites: if Start is at position 80, pick 
            105 (min recommended is 12).""")
            st.write("""
            ***Baseline*** is a region which should be translated, to use as a reference. \n
            ***Target*** is the region you are interested in analysing. You will assess its RTA when using the baseline 
            region activity as reference. \n
            ***Noise Control*** is a region where there should be no translation happening (conceptually akin to a negative 
            control). Use the NOISE SELECTOR tab to identify the right one, unless you have reasons against that.
            """)

    ##### Aggregation Options #####
    with st.expander("1 ~ Select Aggregation Option (optional)"):
        st.write("Here you can select whether to aggregate data (by metadata or overall) or not."
                 "If you want to smooth and aggregate, ALWAYS start from aggregation!")

        valid_agg_options = [None,
                             "Aggregate by Metadata",
                             "Total Aggregation"]

        aggregation_option = st.radio(
            "_",
            label_visibility="collapsed",
            options=valid_agg_options,
            horizontal=True,
            captions=["Default option. No aggregation of the data is applied.",
                      "Aggregates datasets by selected metadata field.",
                      "Aggregates everything, everywhere, all at once."],
        )

        if aggregation_option == "Aggregate by Metadata":
            current_mode.meta_agg_option = st.selectbox(
                "To visualise the diamond plot, select the metadata column to use",
                st.session_state.metadata.columns[:], key="meta_option_aggregation")

        if st.button("Aggregate it!"):
            current_mode.aggregation = aggregation_option
            data_aggregator(current_mode, current_mode.aggregation) # Nb. overwrites all_studies_df (backup still safe)


##### Smoothing Options #####
    with st.expander("2 ~ Select Smoothing Option (optional)"):

        # Failsafe - If no mappability is given, no smoothing is done
        if current_mode.mappability_df is None:
            valid_options = ["No Smoothing", "Criminal Smoothing (Experimental)"]
        else:
            valid_options = ["No Smoothing",
                     "Mappability Smoothing (Proportional)",
                     "Mappability Smoothing (Aggressive)",
                     "Criminal Smoothing (Experimental)"]

        smoothing_option = st.radio(
            "_",
            label_visibility="collapsed",
            options= valid_options,
            horizontal=True,  # makes it look more toggle-like
            captions=["Default option. No smoothing is applied.",
                      "Smoothes the coverage per position by their mappability value (see plot).",
                      "Sets all non-1 mappability positions to 0-counts.",
                      "Per transcript region, finds extremes and smoothes them. Z-score 4 is a good start."],
        )

        if smoothing_option == "Criminal Smoothing (Experimental)":
            z_score_threshold = st.slider("z_score_threshold", 1, 10)
            c1, c2 = st.columns(2)
            with c1:
                st.caption("Abrasive (removes a lot)")
            with c2:
                st.caption("<div style='text-align:right'>Permissive (allows peaks) </div>", unsafe_allow_html=True)

        if st.button("Smooth it!"):
            # Creates backup here?

            current_mode.smoothing = smoothing_option
            if smoothing_option == "Criminal Smoothing (Experimental)":
                current_mode.smoothing_criminal_thresh = z_score_threshold
            else:
                current_mode.smoothing_criminal_thresh = None
            # overwrites all_studies_df (ndr backup still safe) to add smoothing options
            txid_data_smoother(current_mode, current_mode.smoothing)

            after_smooth_plot = plot_ga_with_frames2(current_mode, current_mode.all_studies_df,
                                                     y_column="Global Aggregate smoothed",
                                                     normalize=False)

            streamlit_bokeh(after_smooth_plot, use_container_width=False)




        # Functions that takes in current_mode, meta_agg_option, and does aggregation (stores it too for report).


##### Workflow Recap #####
    st.markdown("#### Recap of your current workflow:")

    col1, col2, col3= st.columns([4, 4, 4])

    with col1: # Ranges Overview
        if current_mode.n_slice is None or current_mode.n_slice[1] == 0:
            bubble("No Ranges Selected", bg_color="#000000", outline_color="#5c5c5c", text_color="#5c5c5c", outline_width=1)
        else:
            bubble("Ranges Selected ", bg_color="#000000", outline_color="#f0aa2d", text_color="#f0aa2d", outline_width=1)

    with col2: # Aggregation Overview
        if current_mode.aggregation is None: # and current_mode.meta_agg_option is None:
            bubble("No Aggregation ", bg_color="#000000", outline_color="#5c5c5c", text_color="#5c5c5c", outline_width=1)
        elif current_mode.aggregation == "Total Aggregation":
            bubble(f"{current_mode.aggregation}", bg_color="#000000",
                   outline_color="#f0aa2d", text_color="#f0aa2d", outline_width=1)
        else:
            bubble(f"{current_mode.aggregation} : {current_mode.meta_agg_option}", bg_color="#000000", outline_color="#f0aa2d", text_color="#f0aa2d", outline_width=1)

    with col3: # Smoothing Overview
        if current_mode.smoothing is None:
            bubble("No Smoothing ", bg_color="#000000", outline_color="#5c5c5c", text_color="#5c5c5c",
                   outline_width=1)
        else:
            bubble(f"{current_mode.smoothing}", bg_color="#000000", outline_color="#f0aa2d", text_color="#f0aa2d", outline_width=1)

    st.write("")

    if st.button("Next"):
        current_mode.txid_ranges = {"x_slice": current_mode.x_slice,
                                   "y_slice": current_mode.y_slice,
                                   "n_slice": current_mode.n_slice}

        st.switch_page("pages/Transcript_Mode_Results.py")