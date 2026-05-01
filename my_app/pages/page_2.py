# page_2.py
import streamlit as st
from app_helpers import *
from module0 import *
from transcriptID_mode_helpers import *
from transcriptID_mode_helpers2 import *
from STP_helpers import *

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

# Renders the info message for this mode
txid_mode_intro()



# Checkpoint: checks for presence of a txid (eg. ENST...) if user has not submitted one, returns error.
txid_start_checkpoint(current_mode, analysis)

# Extract coverage info and transcript sequence. Includes loader.
if current_mode.loading_tx_cov:
    all_ranges, strandedness = txid_data_extraction(current_run = current_mode,
                                                    bigwig_files = analysis.filters.selected_studies,
                                                    gr = st.session_state.gr)

    current_mode.transcript_sequence = get_sequence(fasta, all_ranges, strandedness)

    current_mode.all_studies_df_backup = current_mode.all_studies_df.copy()

# Checks whether the coverage has been extracted, before proceeding
if current_mode.tx_cov_extracted:

    # Format the coverage data as a table for download
    if current_mode.transcript_sequence != "" and analysis.tx_upload == False:
        format_txid_df_download(current_mode)

    # If the sequence is not available and is "File Upload Mode", retrieve sequence from it.
    if current_mode.transcript_sequence is None and analysis.tx_upload == True:
        st.write("You are currently in: Ensembl Transcript ID - File Upload mode")
        try:
            current_mode.transcript_sequence = ''.join(current_mode.all_studies_df['Sequence'])
            current_mode.all_studies_df.drop(columns='Sequence', inplace=True)
            st.info("Retrieved sequence from the uploaded file!")
        except KeyError:
            st.error("The File Uploaded does not come from RiboMARTA ('Sequence' column is missing).")


### Display the data once extracted | Allow for ranges selection ###
if current_mode.tx_cov_extracted == True and current_mode.transcript_sequence != "":

    # # computes global aggregate
    # ga_df_norm = global_aggregate(current_mode.all_studies_df, normalize=True)
    #
    # # Generate plots if they don't exist (lazy initialization via direct checks)
    # if current_mode.tx_ga_plot is None:
    #     # Single transcript plot of global aggregate. No frames or ORF plots
    #     current_mode.tx_ga_plot = plot_global_aggregate_with_zero(ga_df_raw,
    #                                                              current_mode.utr5_length,
    #                                                              current_mode.cds_length,
    #                                                              current_mode.txid)

    if current_mode.triplet_ga_plot is None:
        current_mode.triplet_ga_plot = plot_ga_with_frames2(current_mode,
                                                            current_mode.all_studies_df,
                                                            y_column = "Global Aggregate smoothed",
                                                            normalize=False)

    if current_mode.norm_triplet_ga_plot is None:
        current_mode.norm_triplet_ga_plot = plot_ga_with_frames2(current_mode,
                                                                 current_mode.all_studies_df,
                                                                 y_column = "Global Aggregate smoothed",
                                                                 normalize=True)

    # Expanders to render the plots, if present
    ga_plots_UI_render(current_mode, analysis)


    # Ranges selection for analysis
    # if st.button("Select Ranges"):
    #     current_mode.show_ranges_selectors = True
    with st.expander(" Select Ranges"):

    #if current_mode.show_ranges_selectors:
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
            control). Pick a region with good mappability in the 3'UTR, possibly between stop codons.
            """)

        #st.markdown("<br><br>", unsafe_allow_html=True)

    with st.expander("Select Smoothing Option"):
        # col1, col2, col3 = st.columns([5, 5, 5])
        # with col1:
        #     on = st.toggle("Activate feature")
        #     if st.button("No Smoothing (App Default choice)", key = "no_smoothing"):
        #         current_mode.smoothing_option = None
        #     if st.button("Mappability Smoothing (Proportional)", key = "map_prop_smoothing"):
        #         current_mode.smoothing_option = "Mappability - Proportional"
        #     if st.button("Mappability Smoothing (Aggressive)", key = "map_agg_smoothing"):
        #         current_mode.smoothing_option = "Mappability - Aggressive"

        choice = st.radio(
            "",
            options=["No Smoothing",
                     "Mappability Smoothing (Proportional)",
                     "Mappability Smoothing (Aggressive)",
                     "Criminal Smoothing (Experimental)"],
            horizontal=True,  # makes it look more toggle-like
            captions=["Default option. No smoothing is applied.",
                      "Smoothes the coverage per position by their mappability value (see plot).",
                      "Considers all non-1 mappability positions as 0. Aggressive.",
                      "Per transcript region, finds extreme positions and smoothes them. Z-score 4 is a good start."],
        )


        if choice == "Criminal Smoothing (Experimental)":
            z_score_threshold = st.slider("z_score_threshold", 1, 10)
            c1, c2 = st.columns(2)
            with c1:
                st.caption("Abrasive (removes a lot)")
            with c2:
                st.caption("<div style='text-align:right'>Permissive (allows peaks) </div>", unsafe_allow_html=True)

        if st.button("Restore Data (un-smooth)"):
            if current_mode.all_studies_df.equals(current_mode.all_studies_df_backup):
                st.write("Your dataframe is already consisting of raw counts!")
            else:
                current_mode.all_studies_df = current_mode.all_studies_df_backup.copy()
                st.write("Data has been restored to raw counts!")
        # Add command to backup copy all_studies_df
        # Needs button to reset smoothing! (reassign original all_studies_df


    if st.button("Next"):
        current_mode.txid_ranges = {"x_slice": current_mode.x_slice,
                                   "y_slice": current_mode.y_slice,
                                   "n_slice": current_mode.n_slice}

        current_mode.smoothing_option = choice
        current_mode.smoothing_criminal_thresh = z_score_threshold


        st.switch_page("pages/page_5.py")

