import streamlit as st
from app_helpers import load_config, retrieve_config_path

# Layout
st.set_page_config(layout="wide")

# Load configuration file - default to app internal one if none is given
if "config" not in st.session_state:
    config_path = retrieve_config_path()
    st.session_state["config"] = load_config(config_path)

# Name each page in the navigation bar
home = st.Page("pages/Home.py", title="Home - Organism/Analysis Setup", icon="🏠")
step_0 = st.Page("pages/Metadata_Selector.py", title="Metadata - Dataset Selector", icon="🔍")
step_1_genco = st.Page("pages/Genomic_Mode_Setup.py", title="Genomic Mode - Setup", icon="🧭")
step_2_genco = st.Page("pages/Genomic_Mode_Results.py", title="Genomic Mode - Results", icon="📑")
step_1_txid_mode = st.Page("pages/Transcript_Mode_Setup.py", title="Transcript Mode - Setup", icon="🧭")
step_2_txid_mode = st.Page("pages/Transcript_Mode_Results.py", title="Transcript Mode - Results", icon="📑")
step_3_downloads = st.Page("pages/Results_Export.py", title="Report Download", icon="💾")

# Set navigation bar with my own order of pages
pg = st.navigation([home, step_0, step_1_genco, step_2_genco, step_1_txid_mode, step_2_txid_mode, step_3_downloads])
pg.run()