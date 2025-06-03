# page_2.py
import streamlit as st
from module0 import *
from transcriptID_mode_helpers import *
import os
import pandas as pd
import numpy as np

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


# Required Resources (should be placed into a config.file)
#bigwigs_dir = "/home/gionmattia/Desktop/OAZ3/bed_and_bigwigs/m_musculus/a_shifted_mus_musculus/a_shifted"
#bigwig_files = files = [os.path.abspath(os.path.join(bigwigs_dir, f)) for f in os.listdir(bigwigs_dir) if os.path.isfile(os.path.join(bigwigs_dir, f))]
#gtf_path = "/home/gionmattia/Desktop/OAZ3/annotations/Mus_musculus.GRCm39.113.gtf"


# Sanity check for the transcript ID
txid = st.session_state.get("txid", None)
if txid is None or txid == "":
    col1, col2 = st.columns([1,  3])
    with col1:
        st.image("martas/error_marta.png", width=127)
    with col2:
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.warning(" The 'Transcript ID' field was empty. Please provide a value before proceeding!")

# calls the execution to extract from a single transcript.
else:
    if st.button("Load Transcript"):
        # this sets up the page
        col1, col2, col3 = st.columns([1, 4, 1])
        with col1:
            # Create a placeholder for the image (yes Anmol, I want a picture while results load)
            struggling_image = "martas/struggling_marta.png"
            img_placeholder = st.empty()
            img_placeholder.image(struggling_image, width=100)
        with col2:
            st.markdown("<br><br>", unsafe_allow_html=True)
            with st.spinner("I am retrieving the coverage data..."):
                #st.image("struggling_marta.png", width=100)

                # 1 Retrieve the pyrange object
                try:
                    transcript_ranges = extract_tx_pyranges(st.session_state.gtf_path, st.session_state.txid)

                except KeyError:
                    st.error("I am sorry, but the Transcript ID given was not found in the annotation file."
                             "Please, restart from the Home page.")
                    img_placeholder.image("martas/ko_marta.png", width = 100)
                    st.stop()
                    # Is there a way to stop the spinner here?

                # 2  Retrieve features (CDS, UTRs) lengths
                all_ranges, _ = regions_dataframe_formatter(transcript_ranges)
                st.session_state.utr5_length, st.session_state.cds_length, st.session_state.utr3_length = extract_features_length(all_ranges)
                # 3 Extract coverage from given bigwigs
                st.session_state.all_studies_df = extract_from_multiple_bigwigs(st.session_state.bigwig_files,
                                                                                transcript_ranges,
                                                                                n_jobs=2)

                # Resets the picture from above
                img_placeholder.empty()

        # Allow to save the raw data as a .csv
        csv = st.session_state.all_studies_df.to_csv(index=True).encode('utf-8')

        st.download_button(
            label="Download data as CSV",
            data=csv,
            file_name= f'{st.session_state.txid}_raw_counts.csv',
            mime='text/csv'
        )

    if 'picture' not in st.session_state:
        # computes global aggregate
        if st.session_state.get("all_studies_df") is None:
            st.stop()
        ga_df = global_aggregate(st.session_state.all_studies_df)
        # creates single transcript plot
        st.session_state.picture = plot_global_aggregate_with_zero(ga_df,
                                                                   st.session_state.utr5_length,
                                                                   st.session_state.cds_length,
                                                                   st.session_state.txid)
    # if picture already exists, renders it
    if st.session_state.get("picture") is not None:
        st.plotly_chart(st.session_state['picture'], use_container_width=True)

    # It would be best if the plotly also had an ORF plot beneath...maybe later on.

    if "show_inputs" not in st.session_state:
        st.session_state.show_inputs = False
    if st.button("Select Ranges"):
        st.session_state.show_inputs = True

    if st.session_state.show_inputs:
        col1, col2, col3 = st.columns([5, 5, 5]) # originally [5, 1, 5]
        full_mRNA_length = st.session_state.utr5_length + st.session_state.cds_length + st.session_state.utr3_length

        with col1:
            st.write("Select First Region")
            st.session_state.x_slice = ranges_selector(full_mRNA_length, "x")

        with col2:
            st.write("Select Second Region")
            st.session_state.y_slice = ranges_selector(full_mRNA_length, "y")

        with col3:
            st.write("Select Noise Region")
            st.session_state.n_slice = ranges_selector(full_mRNA_length, "n")

        if st.button("Next"):
            # NEED to add an internal check here
            st.session_state.coordinates = {"x_slice": st.session_state.x_slice,
                                            "y_slice": st.session_state.y_slice,
                                            "n_slice": st.session_state.n_slice}
            st.switch_page("pages/page_4.py")
