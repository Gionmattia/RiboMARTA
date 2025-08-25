# page_2.py
import streamlit as st
from app_helpers import *
from module0 import *
from transcriptID_mode_helpers import *
from transcriptID_mode_helpers2 import *

#First checkpoint. Stops rendering if the user has not initialised any analysis in the Home page
quirky_id_checkpoint()

# Retrieve the AnalysysStatus class corresponding to the quirky_id, from the dictionary "runs"
current_run = st.session_state.runs[st.session_state.quirky_id]
fasta = st.session_state.fasta
st.write(f"######  :gray[ Analysis ID: {st.session_state.quirky_id}]")
st.title("Transcript ID Mode")

# Renders the info message for this mode
txid_mode_intro()

# Checkpoint: checks for presence of a txid (eg. ENST...) if user has not submitted one, returns error.
txid_start_checkpoint(current_run)

# Extract coverage info and transcript sequence. Includes loader.
if current_run.loading_tx_cov:
    all_ranges, strandedness = txid_data_extraction(current_run, st.session_state.gr)
    current_run.transcript_sequence = get_sequence(fasta, all_ranges, strandedness)


# Checks whether the coverage has been extracted, before proceeding
if current_run.tx_cov_extracted:

    # Format the coverage data as a table for download
    if current_run.transcript_sequence != "" and current_run.tx_upload == False:
        format_txid_df_download(current_run)

    # If the sequence is not available and is "File Upload Mode", retrieve sequence from it.
    if current_run.transcript_sequence == "" and current_run.tx_upload == True:
        st.write("You are currently in: Ensembl Transcript ID - File Upload mode")
        try:
            current_run.transcript_sequence = ''.join(current_run.all_studies_df['Sequence'])
            current_run.all_studies_df.drop(columns='Sequence', inplace=True)
        except KeyError:
            st.error("The File Uploaded does not come from RiboMARTA ('Sequence' column is missing).")


# For the next chunk, we need coverage AND transcript sequence to be present
if current_run.tx_cov_extracted == True and current_run.transcript_sequence != "":

    # computes global aggregate
    ga_df = global_aggregate(current_run.all_studies_df) # SRR[...i] | "Global Aggregate"
    # creates single transcript plot for it
    current_run.tx_ga_plot = plot_global_aggregate_with_zero(ga_df,
                                                             current_run.utr5_length,
                                                             current_run.cds_length,
                                                             current_run.txid)

    # creates single transcript plot with frame info
    if current_run.triplet_ga_plot is None:
        make_triplet_ga_plot(current_run, ga_df)

    # Expanders to render the plots, if present
    ga_plots_UI_render(current_run)

    if st.button("Select Ranges"):
        current_run.show_ranges_selectors = True

    if current_run.show_ranges_selectors:
        col1, col2, col3 = st.columns([5, 5, 5])
        full_mRNA_length = current_run.utr5_length + current_run.cds_length + current_run.utr3_length

        with col1:
            st.write("Select First Region")
            current_run.x_slice = ranges_selector(full_mRNA_length, "x")

        with col2:
            st.write("Select Second Region")
            current_run.y_slice = ranges_selector(full_mRNA_length, "y")

        with col3:
            st.write("Select Noise Region")
            current_run.n_slice = ranges_selector(full_mRNA_length, "n")

        st.markdown("<br><br>", unsafe_allow_html=True)

        # Renders info on the Frame Aware mode.
    frame_aware_analysis_intro()

    col1, col2 = st.columns([3, 3])
    with col1:
        current_run.frame_aware = st.toggle("Perform Frame-Aware Analysis", key = "frame_aware_selector")
    with col2:
        # Add a Marta with a baloon specifying the analysis selected!
        if current_run.frame_aware:
            st.write("Frame-Aware Analysis Activated!")
        elif not current_run.frame_aware:
            st.write("Classic Analysis Selected!")

    st.markdown("<br><br>", unsafe_allow_html=True)

    if st.button("Next"):

        if not current_run.frame_aware:
            current_run.txid_ranges = {"x_slice": current_run.x_slice,
                                       "y_slice": current_run.y_slice,
                                       "n_slice": current_run.n_slice}
            st.switch_page("pages/page_5.py")

        elif current_run.frame_aware:
            current_run.txid_ranges = {"x_slice": current_run.x_slice,
                                       "y_slice": current_run.y_slice,
                                       "n_slice": current_run.n_slice}
            st.switch_page("pages/page_6.py")