# page_2.py
import streamlit as st
from smoothing_helpers import *
from app_helpers import *
from module0 import *
from transcriptID_mode_helpers import *
from transcriptID_mode_helpers2 import *
from STP_helpers import *
from streamlit_bokeh import streamlit_bokeh

import streamlit as st

st.markdown("""
<style>
    /* 1. The Base Button Style (Targeting 1st and 2nd) */
    div[data-testid="stButton"]:nth-of-type(1) > button,
    div[data-testid="stButton"]:nth-of-type(2) > button,
    div[data-testid="stDownloadButton"] > button {
        background-color: #000000;
        color: #ffffff;
        border: 2px solid #555555;
        border-radius: 10px;
        padding: 0.5em 1em;
        transition: all 0.3s ease;
    }

    /* 2. The Hover State (Targeting 1st and 2nd) */
    div[data-testid="stButton"]:nth-of-type(1) > button:hover,
    div[data-testid="stButton"]:nth-of-type(2) > button:hover,
    div[data-testid="stDownloadButton"] > button:hover  {
        border-color: #FFA500;
        color: #FFA500;
        background-color: #1a1a1a;
    }

    /* 3. The Clicked/Active State (Targeting 1st and 2nd) */
    div[data-testid="stButton"]:nth-of-type(1) > button:active,
    div[data-testid="stButton"]:nth-of-type(2) > button:active,
    div[data-testid="stDownloadButton"] > button:active {
        background-color: #FFA500;
        color: #ffffff;
        border-color: #FFA500;
    }

    /* 4. Focus State (Targeting 1st and 2nd) */
    div[data-testid="stButton"]:nth-of-type(1) > button:focus:not(:active),
    div[data-testid="stButton"]:nth-of-type(2) > button:focus:not(:active),
    div[data-testid="stDownloadButton"] button:focus:not(:active) >
     {
        border-color: #FFA500;
        color: #FFA500;
        box-shadow: none;
    }

        

/* ===== EXPANDERS ===== */
/* New versions */
div[data-testid="stExpander"] {
    background-color: #2b2b2b !important;
    border-radius: 8px;
    border: 1px solid #444;
}

div[data-testid="stExpander"] > details {
    background-color: #f0aa2d !important;
}

div[data-testid="stExpander"] > details > summary,
div[data-testid="stExpander"] > details > summary p,
div[data-testid="stExpander"] > details > summary span {
    background-color: #f0aa2d !important;
    color: #000000 !important;
    padding: 0.17rem 0.5rem;
    font-size: 1.15rem !important;
}

div[data-testid="stExpander"] > details[open] > summary {
    background-color: #f0aa2d !important;
}

/* Content inside expander */
div[data-testid="stExpanderDetails"] {
    background-color: #2b2b2b !important;
    padding: 1rem;
}

/* ===== TABS ===== */
button[data-baseweb="tab"] {
    background-color: #f0aa2d !important;
    color: #000000 !important;
    padding: 0.75rem 1.5rem !important;
    font-size: 1.15rem !important;
    font-weight: 600 !important;
}
button[data-baseweb="tab"] p {
    font-size: 1.1rem !important;
}

button[data-baseweb="tab"][aria-selected="true"] {
    background-color: #f0aa2d !important;
    color: #000000 !important;
}

div[data-baseweb="tab-highlight"] {
    background-color: #f0c694 !important;
    height: 0.45rem !important;
}

div[data-baseweb="tab-panel"] {
    background-color: #2b2b2b !important;
    padding: 1.5rem;
    border-radius: 0 0 8px 8px;
}

div[data-baseweb="tab-list"] {
    background-color: transparent !important;
    gap: 4px;
}
/* Unselected radio - default state (outer ring) */
div[data-testid="stRadio"] label[data-baseweb="radio"] > div:first-child {
    border: 2px solid #888888 !important;
    background-color: transparent !important;
    box-sizing: border-box !important;
}
/* Unselected radio - hover state (outer ring) */
div[data-testid="stRadio"] label[data-baseweb="radio"]:hover:not(:has(input:checked)) > div:first-child {
    border: 2px solid #d0a454 !important;
    background-color: rgba(30, 144, 255, 0.1) !important;
    box-sizing: border-box !important;
}
/* Outer ring of the selected radio */
div[data-testid="stRadio"] label[data-baseweb="radio"]:has(input:checked) > div:first-child {
    border: 2px solid #d0a454 !important;
    background-color: transparent !important;
    box-sizing: border-box !important;
}
/* Inner dot of the selected radio */
div[data-testid="stRadio"] label[data-baseweb="radio"]:has(input:checked) > div:first-child > div {
    border: 2px solid #d0a454 !important;
    background-color: #d0a454 !important;
    box-sizing: border-box !important;
}
</style>
""", unsafe_allow_html=True)


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
# Checkpoint: checks for presence of a txid (eg. ENST...) if user has not submitted one, returns error.
txid_checkpoint(current_mode)
# Checkpoint - allows to extract coverage only if we are not in "Upload Mode"
if analysis.tx_upload == False:
    if st.button("Load Transcript"):  # If not in tx upload, loads
        current_mode.loading_tx_cov = True

## All the code to extract info - needs some order ##
# Extract coverage info, transcript sequence, and mappability. Includes loader.
if current_mode.loading_tx_cov:
    # Extracts mappability info, if available.
    if analysis.mappability_path is not None:
        mappability_df = txid_map_extraction(current_mode.txid, analysis.mappability_path, st.session_state.gr)
        mappability_df["X"] = mappability_df.index
        other_col = [c for c in mappability_df.columns if c not in ['Position', 'X']][0]
        current_mode.mappability_df = mappability_df.rename(columns={other_col: 'Mappability'})

    current_mode.aggregation_option = ""
    current_mode.smoothing_option = False
    all_ranges, strandedness = txid_data_extraction(current_run = current_mode,
                                                    bigwig_files = analysis.filters.selected_studies,
                                                    gr = st.session_state.gr)

    current_mode.transcript_sequence = get_sequence(fasta, all_ranges, strandedness)

    # Creates back-ups for the smoothing options and aggregation options.
    current_mode.all_studies_df_backup = current_mode.all_studies_df.copy()
    current_mode.all_studies_df_pre_aggregation_backup = current_mode.all_studies_df.copy()


# Checks whether the coverage has been extracted, before proceeding
if current_mode.tx_cov_extracted:
    # If data has been extracted, allows for download in .csv format
    if analysis.tx_upload == False:
        format_txid_df_download(current_mode)
    # If the sequence is not available and is "File Upload Mode", retrieve sequence from it.
    elif analysis.tx_upload == True:
        st.markdown("### You are currently in: Ensembl Transcript ID - File Upload mode")
else:
    st.stop()


### Display the data once extracted | Allow for ranges selection ###
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


    with st.container(border = True): #"Select Ranges"
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


##### Smoothing Options #####
    with st.expander("Option A ~ Select Smoothing Option (optional)"):

        # Failsafe - If no mappability is given, no smoothing is done
        if current_mode.mappability_df is None:
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
            # Users cannot AGGREGATE and SMOOTH at the same time. One option negates the other.
            current_mode.aggregation_option = None
            current_mode.meta_agg_option = None
            # Undoes any aggregation does beforehand
            data_aggregator(current_mode, current_mode.aggregation_option)

            current_mode.smoothing_option = choice
            if choice == "Criminal Smoothing (Experimental)":
                current_mode.smoothing_criminal_thresh = z_score_threshold
            else:
                current_mode.smoothing_criminal_thresh = None
            # overwrites all_studies_df (ndr backup still safe) to add smoothing options
            txid_data_smoother(current_mode, current_mode.smoothing_option)

            after_smooth_plot = plot_ga_with_frames2(current_mode, current_mode.all_studies_df,
                                                     y_column="Global Aggregate smoothed",
                                                     normalize=False)

            streamlit_bokeh(after_smooth_plot, use_container_width=False)


##### Aggregation Options #####
    with st.expander("Option B ~ Select Aggregation Option (optional)"):
        st.write("Here you can select whether to aggregate data (by metadata or overall) or not. Resetting this options"
                 "resets also any smoothing option previously selected!")

        valid_agg_options = [None,
                         "Aggregate by metadata field",
                         "Total Aggregation"]

        aggregation_option = st.radio(
            "_",
            label_visibility="collapsed",
            options=valid_agg_options,
            horizontal=True,  # makes it look more toggle-like
            captions=["Default option. No aggregation of the data is applied.",
                      "Aggregates datasets by selected metadata field.",
                      "Aggregates everything, everywhere, all at once."],
        )

        if aggregation_option == "Aggregate by metadata field":
            current_mode.meta_agg_option = st.selectbox(
                "To visualise the diamond plot, select the metadata column to use",
                st.session_state.metadata.columns[:], key =  "meta_option_aggregation")

        if st.button("Aggregate it!"):
            # Users cannot AGGREGATE and SMOOTH at the same time. One option negates the other.
            current_mode.smoothing_option = False
            # Undoes any smoothing
            txid_data_smoother(current_mode, current_mode.smoothing_option)

            current_mode.aggregation_option = aggregation_option
            # Nb. overwrites all_studies_df (ndr backup still safe)
            data_aggregator(current_mode, current_mode.aggregation_option)


        # Functions that takes in current_mode, meta_agg_option, and does aggregation (stores it too for report).


##### Workflow Recap #####
    st.markdown("#### Recap of your current workflow:")

    col1, col2, col3= st.columns([4, 4, 4])

    with col1:
        if current_mode.n_slice is None or current_mode.n_slice[1] == 0:
            bubble("No Ranges Selected", bg_color="#000000", outline_color="#5c5c5c", text_color="#5c5c5c", outline_width=1)
        else:
            bubble("Ranges Selected ", bg_color="#000000", outline_color="#f0aa2d", text_color="#f0aa2d", outline_width=1)
    with col2:
        if current_mode.smoothing_option is False:
            bubble("No Smoothing ", bg_color="#000000", outline_color="#5c5c5c", text_color="#5c5c5c",
                   outline_width=1)
        else:
            bubble(f"{current_mode.smoothing_option}", bg_color="#000000", outline_color="#f0aa2d", text_color="#f0aa2d", outline_width=1)
    with col3:
        if current_mode.aggregation_option is None and current_mode.meta_agg_option is None:
            bubble("No Aggregation ", bg_color="#000000", outline_color="#5c5c5c", text_color="#5c5c5c",
                   outline_width=1)
        elif current_mode.aggregation_option == "Total Aggregation":
            bubble(f"{current_mode.aggregation_option}", bg_color="#000000",
                   outline_color="#f0aa2d", text_color="#f0aa2d", outline_width=1)
        else:
            bubble(f"{current_mode.aggregation_option} : {current_mode.meta_agg_option}", bg_color="#000000", outline_color="#f0aa2d", text_color="#f0aa2d", outline_width=1)

    #: str = None
    st.write("")

    if st.button("Next"):
        current_mode.txid_ranges = {"x_slice": current_mode.x_slice,
                                   "y_slice": current_mode.y_slice,
                                   "n_slice": current_mode.n_slice}

        st.switch_page("pages/page_5.py")

