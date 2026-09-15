# Functions to validate the transcript upload should go down here, ideally.

import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin
import streamlit as st
from home_helpers import *

#  Needs to be implemented as an if/else for that specific hosting service...
# BUT IT COULD WORK
def get_bigwig_urls(base_url):
    urls = []
    try:
        # Connect to server
        response = requests.get(base_url)
        response.raise_for_status()

        # Read the HTML directory listing that RangeHTTPServer creates
        soup = BeautifulSoup(response.text, 'html.parser')

        # Find every link on that webpage
        for link in soup.find_all('a'):
            href = link.get('href')

            # Only keep the BigWig files
            if href and href.endswith('.bw'):
                # urljoin to safely combine
                full_url = urljoin(base_url, href)
                urls.append(full_url)

    except requests.exceptions.RequestException as e:
        st.error(f"I am sorry but I could not connect to the GWIPS server: {e}")

    return urls



# Organism Loader - Provides interface and code to load organism.
# This is the version that checks whether an url or a path was given to load the BigWigs from.
# Substitute to orgasnism_loader_interfact in Home.py
def organism_loader_interface2():
    """
    Displays the Streamlit UI for selecting and loading an organism's data.

    This function creates a select box and a button. When the user clicks the
    button, it loads all necessary configuration and data files (metadata, paths, etc.)
    into st.session_state.

    Side Effects:
        - Populates st.session_state with data like 'gtf_path', 'metadata', 'fasta', etc.
        - Sets st.session_state.config_loaded = True upon successful loading.
    """

    # Load configuration file
    config = load_config("config.yaml")
    # Load path to annotation file, metadata, and the bigwigs in the app
    organisms_available = list(config.keys())

    # The selectbox determines which organism will be loaded
    selected_key = st.selectbox("Select from list", organisms_available)

    if st.button(f"Load {selected_key}"):
        # 1. Update session state with the selected organism and its specific config
        st.session_state.selected_organism = selected_key
        organism_config = config[st.session_state.selected_organism]

        # 2. Load all paths and data into session state
        with st.spinner(f"Loading data for {st.session_state.selected_organism}..."):
            st.session_state.gtf_path = organism_config["gtf_path"]
            st.session_state.metadata = pd.read_csv(organism_config["metadata_path"])
            st.session_state.mappability_path = organism_config.get("mappability_path")

            bigwigs_dir = organism_config["bigwigs_dir"]
            st.session_state.fasta_path = organism_config["fasta_path"]
            # Parse urls if given https
            if "http" in bigwigs_dir:
                st.session_state.bigwig_files = get_bigwig_urls(bigwigs_dir)

                urls = [
                    f.replace(".bw", "") for f in st.session_state.bigwig_files
                ]

                st.session_state.full_dataset_list = [p.split("/")[-1] for p in urls]

            else:
                st.session_state.bigwig_files = [
                    os.path.abspath(os.path.join(bigwigs_dir, f))
                    for f in os.listdir(bigwigs_dir)
                    if os.path.isfile(os.path.join(bigwigs_dir, f))
                ]

                st.session_state.full_dataset_list = [
                    f.replace(".bw", "") for f in os.listdir(bigwigs_dir)
                    if os.path.isfile(os.path.join(bigwigs_dir, f))
                ]



            # These functions can be time-consuming
            st.session_state.fasta = read_fasta(st.session_state.fasta_path)
            st.session_state.gr = read_gtf_annotation(st.session_state.gtf_path)

        # 3. Set a flag indicating that the configuration is loaded
        st.session_state.config_loaded = True
        st.success(f"Loaded config for {st.session_state.selected_organism}!")