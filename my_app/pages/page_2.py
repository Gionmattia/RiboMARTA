# page_2.py
import streamlit as st
from app_helpers import *
from module0 import *
from transcriptID_mode_helpers import *
import time # to time process execution
from bokeh.models import FactorRange
from bokeh.embed import file_html
from bokeh.resources import CDN
from bokeh.layouts import column

# We are retrieving the AnalysysStatus class corresponding to the quirky_id, from the dictionary "runs"
current_run = st.session_state.runs[st.session_state.quirky_id]
st.write(f"######  :gray[ Analysis ID: {st.session_state.quirky_id}]")
st.title("Transcript ID Mode")

col1, col2 = st.columns([4,1])
with col1:
    st.write("Welcome to the Transcript ID Analysis page. This analysis is broken down in two stages:")
    st.write("""1) By clicking on "Load Transcript", you will be shown a plot with the aggregate Ribo-Seq footprints 
    across all the available studies.""")
    st.write("Since this first step can take time, you will be allowed to save the raw data if you wish so. Click on 'Dowload"
             "Data as CSV'")
    st.write(
        """2) Once the dynamic plot is visualised, you can click on 'Select Ranges' to specify the areas of the transcript
         you want to compare in the analysis.""")
with col2:
    st.image("martas/explaining_marta.png", width = 180)

# Sanity check for the transcript ID
txid = current_run.txid
if txid == "":
    col1, col2 = st.columns([1,  3])
    with col1:
        st.image("martas/error_marta.png", width=127)
    with col2:
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.warning(" The 'Transcript ID' field was empty. Please provide a value before proceeding!")

# Allows to extrat the tx coverage
else:
    if st.button("Load Transcript"):
        current_run.loading_tx_cov = True

if current_run.loading_tx_cov == True:
    # this sets up the page
    col1, col2, col3 = st.columns([1, 4, 1])
    with col1:
        # Create a placeholder for the image (yes Anmol, I want a picture while results load)
        struggling_image = "martas/struggling_marta.png"
        img_placeholder = st.empty()
        img_placeholder.image(struggling_image, width=100)
    with col2:
        # st.markdown("<br><br>", unsafe_allow_html=True)
        with st.spinner("I am retrieving the coverage data..."):
            try:
                transcript_ranges = extract_tx_pyranges(st.session_state.gr, current_run.txid)
                fasta = st.session_state.fasta

            except KeyError:
                st.error("I am sorry, but the Transcript ID given was not found in the annotation file."
                             "Please, restart from the Home page.")
                img_placeholder.image("martas/ko_marta.png", width = 100)
                st.stop()

            # 2  Retrieve transcript features (CDS, UTRs) ranges and strandedness
            all_ranges, strandedness = regions_dataframe_formatter(transcript_ranges)
            # 3a  Retrieve features (CDS, UTRs) lengths
            current_run.utr5_length, current_run.cds_length, current_run.utr3_length = extract_features_length(all_ranges)
            # 3b Extract coverage from given bigwigs
            current_run.all_studies_df = extract_from_multiple_bigwigs(st.session_state.bigwig_files,
                                                                            all_ranges,
                                                                            strandedness,
                                                                            n_jobs=1)
                # Resets the picture from above
            img_placeholder.empty()
            current_run.loading_tx_cov = False
            current_run.tx_cov_extracted = True

# Checks whether the coverage has been extracted, before offering other UI options
if current_run.tx_cov_extracted == True:

    # Allow to save the raw data as a .csv
    csv = current_run.all_studies_df.to_csv(index=True).encode('utf-8')
    st.download_button(
        label="Download data as CSV",
        data=csv,
        file_name=f'{current_run.txid}_raw_counts.csv',
        mime='text/csv'
    )

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
    if current_run.transcript_sequence == "":
        current_run.transcript_sequence = get_sequence(fasta, all_ranges, strandedness) # if it works, change with st.session_state.fasta

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
            st.write("Below you can see a Global aggregate for your data, normalised by transcript counts."
                     "This plot is a simplified visualisation. Easy on the eyes, but less informative.")
            if current_run.tx_ga_plot is not None:
                st.plotly_chart(current_run.tx_ga_plot, use_container_width=True)
        with tab2:
            st.write("Below you can see data for your transcript coloured by frame."
                     "To help the detection of patterns, data has been smoothed by a window of size k."
                     "This visualisation is less intuitive, but more informative.")
            if current_run.triplet_ga_plot is not None:
                # Saving as html, since blokeh v.3+ is not directly compatible with streamlit...
                # Patchy, maybe I can come up with something better... but does it job for now.
                html = file_html(current_run.triplet_ga_plot, CDN, "Transcript_plot_with_ORF")
                st.components.v1.html(html, height=700, width=1600)

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

        if st.button("Next"):
            current_run.txid_ranges = {"x_slice": current_run.x_slice,
                                            "y_slice": current_run.y_slice,
                                            "n_slice": current_run.n_slice}
            st.switch_page("pages/page_5.py")
