# Home.py

import streamlit as st
import pandas as pd
from home_helpers import *
import os


# st.session_state.gtf_path = config["gtf_path"]
# st.session_state.metadata = pd.read_csv(config["metadata_path"])
# bigwigs_dir = config["bigwigs_dir"]
# st.session_state.fasta_path = config["fasta_path"]
# st.session_state.bigwig_files = files = [os.path.abspath(os.path.join(bigwigs_dir, f)) for f in os.listdir(bigwigs_dir) if os.path.isfile(os.path.join(bigwigs_dir, f))]


st.title("Welcome to RiboMARTA")
col1, col2 = st.columns([1,4])
with col1:
    st.image("martas/marta.png", width = 200)
with col2:
    st.write("""Welcome to the first draft of RiboMARTA. My name is Marta, and I am here to guide you
    through the steps of the analysis.""")
    
    st.write("""If Pasha likes me, I will be allowed to live.""")
    
    st.write("""You are currently in the Home page of the analysis. In here you can decide whether to extract regions using
    genomic coordinates, an ensembl transcript ID, or use raw data from a previous run fo RiboMARTA. Once you are settled, press
    'Next'.""")


## THIS COULD BE WRAPPED UP IN A SINGLE FUNCTION WITHIN HOME_HELPERS!!
# Load configuration file
config = load_config("config.yaml")
# Load path to annotation file, metadata, and the bigwigs in the app
organisms_available = list(config.keys())

st.write("xx Configure Organism for the analysis:")
st.session_state.selected_organism = st.selectbox("Choose an organism:", organisms_available)
if st.button("Load this organism config"):
    config = config[st.session_state.selected_organism]
    st.session_state.gtf_path = config["gtf_path"]
    st.session_state.metadata = pd.read_csv(config["metadata_path"])
    bigwigs_dir = config["bigwigs_dir"]
    st.session_state.fasta_path = config["fasta_path"]
    st.session_state.bigwig_files = [
        os.path.abspath(os.path.join(bigwigs_dir, f))
        for f in os.listdir(bigwigs_dir)
        if os.path.isfile(os.path.join(bigwigs_dir, f))
    ]

    st.session_state.fasta = read_fasta(st.session_state.fasta_path)
    st.session_state.gr = read_gtf_annotation(st.session_state.gtf_path)
    st.session_state.config_loaded = True
    st.success(f"Loaded config for {st.session_state.selected_organism}!")

# Sets up possible choices for analysis
df = pd.DataFrame({
  'Analysis_options': ["None Selected", "Genomic Coordinates", "Ensembl Transcript ID", "File Upload"]
})

option = st.selectbox(
    'Select the Analysis Mode',
     df['Analysis_options'])

# Renders the different options
if option == 'Genomic Coordinates':
    genomic_coordinates_selection()

elif option == 'Ensembl Transcript ID':
    ensembl_txid_selection()

elif option ==  "File Upload":
    upload_selection()

elif option ==  "None Selected":
    st.write("Please select the analysis mode you desire.")
    # session_state clean-up (in this case, removes everything)
    keys_to_reset = ["my_ranges", "picture", "txid", "genomic_slices", "all_studies_df", "utr5_length",
                             "cds_length", "utr3_length"]
    reset_session_state(keys_to_reset)



