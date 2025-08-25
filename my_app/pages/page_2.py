# page_2.py
import streamlit as st
from app_helpers import *
from module0 import *
from transcriptID_mode_helpers import *
from transcriptID_mode_helpers2 import *
from bokeh.models import FactorRange
from bokeh.embed import file_html
from bokeh.resources import CDN
from bokeh.layouts import column

#First checkpoint. Stops rendering if the user has not initialised any analysis in the Home page
quirky_id_checkpoint()

# Retrieve the AnalysysStatus class corresponding to the quirky_id, from the dictionary "runs"
current_run = st.session_state.runs[st.session_state.quirky_id]
fasta = st.session_state.fasta
st.write(f"######  :gray[ Analysis ID: {st.session_state.quirky_id}]")
st.title("Transcript ID Mode")

# renders the info message for this mode
txid_mode_intro()

# Checkpoint: checks for presence of a txid (eg. ENST...) if user has not submitted one, returns error.
txid_start_checkpoint(current_run)

# Extract coverage info. Includes loader.
if current_run.loading_tx_cov == True:
    all_ranges, strandedness = txid_data_extraction(current_run, st.session_state.gr)


# Checks whether the coverage has been extracted, before proceeding
if current_run.tx_cov_extracted == True:

    # If the txid sequence is not available, retrieve it.
    if current_run.transcript_sequence == "" and current_run.tx_upload == False:
        current_run.transcript_sequence = get_sequence(fasta, all_ranges, strandedness)

    # If the sequence is retrieved, add it to the dataframe before download.
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
    # Calculates frame and smoothes data (smoothing ny window size k)
    smoothed_ga_df, k = smoother(ga_df)
    # Defines x-axis as FactorRange object (allows for multiple plots to share X)
    positions = smoothed_ga_df['X'].astype(str).tolist()
    x_range = FactorRange(*positions)

    # --- Plotting Functions --- #

    # creates single transcript plot
    current_run.tx_ga_plot = plot_global_aggregate_with_zero(ga_df,
                                                             current_run.utr5_length,
                                                             current_run.cds_length,
                                                             current_run.txid)

    # creates single transcript plot with frame info
    if current_run.triplet_ga_plot is None:
        picture2_transcript_plot = subplot_single_transcript_plot(smoothed_ga_df,
                                                                  'Global_Aggregate_smoothed',
                                                                  x_range,
                                                                  k,
                                                                  start = current_run.utr5_length, # NEED TO FIX THIS
                                                                  stop = current_run.cds_length + current_run.utr5_length,
                                                                  txid = current_run.txid)

        picture2_transcript_plot.height = 400  # or whatever height you need
        picture2_transcript_plot.sizing_mode = 'stretch_width'

        # Code to retrieve sequence and codon positions
        orf_info = compute_codons_by_frame(current_run.transcript_sequence)
        # Code to create the ORF plots and the combined plot.
        orfplot_f1 = ORF_plotter(orf_info, 1, x_range)
        orfplot_f1.height = 75
        orfplot_f1.sizing_mode = 'stretch_width'
        orfplot_f2 = ORF_plotter(orf_info, 2, x_range)
        orfplot_f2.height = 75
        orfplot_f2.sizing_mode = 'stretch_width'
        orfplot_f3 = ORF_plotter(orf_info, 3, x_range)
        orfplot_f3.height = 75
        orfplot_f3.sizing_mode = 'stretch_width'
        current_run.triplet_ga_plot = column(picture2_transcript_plot, orfplot_f1, orfplot_f2, orfplot_f3, sizing_mode='stretch_width')

    # # if picture already exists, renders it
    with st.expander("Global Aggregate Plots"):
        tab1, tab2 = st.tabs(["Simple Plot", "Periodicity Plot"])
        with tab1:
            st.write("Below you can see a Global aggregate for your data (the total sum of raw counts from each study)."
                     "This plot is a simplified visualisation. Easy on the eyes, but less informative for what concernes"
                     "the possible reading frames.")
            if current_run.tx_ga_plot is not None:
                st.plotly_chart(current_run.tx_ga_plot, use_container_width=True)

        with tab2:
            st.write("Below you can see a Global Aggregate for the (normalised) data, coloured by frame."
                     "To help the detection of patterns, data has been smoothed by a window of size k."
                     "This visualisation is less intuitive (since the y-axis now displays a normalised value),"
                     " but more informative for what concerns the possible reading frames translated.")
            if current_run.triplet_ga_plot is not None:
                # Saving as html, since blokeh v.3+ is not directly compatible with streamlit...
                # Patchy, maybe I can come up with something better... but does it job for now.
                html = file_html(current_run.triplet_ga_plot, CDN, "Transcript_plot_with_ORF")
                st.components.v1.html(html, height=650, width=1600)
            else:
                st.write("triplet_GA_missing")

    if "show_ranges_selectors" not in st.session_state:
        st.session_state.show_ranges_selectors = False
    if st.button("Select Ranges"):
        st.session_state.show_ranges_selectors = True

    if st.session_state.show_ranges_selectors:
        col1, col2, col3 = st.columns([5, 5, 5]) # originally [5, 1, 5]
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

        st.write("### About Frame Awareness")
        st.write("The analysis uses coverage data between regions coming from all three different reading frames."
                 "While this allows to have stronger signal, it is less precise and cannot allow for the discrimination"
                 "of overlapping or nested ORFs.")
        st.write("The Frame-Aware analysis allows to compute the RTA using signal coming from each reading frame "
                 "individually. While this is mode is more precise, it has the downside not all individual studies will"
                 " have an in-frame signal strong enough to be detected as significant.")
        st.write("")

        col1, col2 = st.columns([3, 3])
        with col1:
            current_run.frame_aware = st.toggle("Perform Frame-Aware Analysis", key = "frame_aware_selector")
        with col2:
            # Add a Marta with a baloon specifying the analysis selected!
            if current_run.frame_aware == True:
                st.write("Frame-Aware Analysis Activated!")
            elif current_run.frame_aware == False:
                st.write("Classic Analysis Selected!")

        st.markdown("<br><br>", unsafe_allow_html=True)

        if st.button("Next"):
            if current_run.frame_aware == False:
                current_run.txid_ranges = {"x_slice": current_run.x_slice,
                                           "y_slice": current_run.y_slice,
                                           "n_slice": current_run.n_slice}
                st.switch_page("pages/page_5.py")

            elif current_run.frame_aware == True:
                current_run.txid_ranges = {"x_slice": current_run.x_slice,
                                           "y_slice": current_run.y_slice,
                                           "n_slice": current_run.n_slice}
                st.switch_page("pages/page_6.py")

# ALSO, REWORK PICTURES SO THEY ARE SAVED AND NOT RECOMPUTED!!