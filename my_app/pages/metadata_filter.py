import streamlit as st
import pandas as pd
import os

analysis = st.session_state.analysis_runs[st.session_state.quirky_id]
st.write(f"######  :gray[ Analysis ID: {analysis.quirky_id}]")
st.title("Dataset Selection")
st.write("Here you can select which datasets to process, based on the available metadata. "
         "The current interface allows to choose from any metadata information available on the RDP,"
         " good suggestions to start from are TISSUE, CONDITION, or CELL_LINE")

# Move the dataset list to the current analysis run
analysis.full_dataset_list = st.session_state.full_dataset_list
# Retrieve absolute paths for all bigwigs
bigwigs = st.session_state.bigwig_files
# Subset the metadata based on the availability of bigwigs.
restricted_meta = st.session_state.metadata[st.session_state.metadata["Run"].isin(analysis.full_dataset_list)]
# Create selector based on available metadata
available_meta =  restricted_meta.columns
selected_columns = st.multiselect(f"Select the metadata fields to filter by", available_meta)

# Build filters dictionary dynamically
filters = {}
for col in selected_columns:
    unique_vals = restricted_meta[col].unique()
    selected = st.multiselect(f"Select values for {col}", unique_vals)
    if selected:
        filters[col] = selected

# Apply all filtering criteria
mask = pd.Series([True] * len(restricted_meta), index=restricted_meta.index)
for col, values in filters.items():
    mask &= restricted_meta[col].isin(values)

filtered_df = restricted_meta[mask]

# Message to user
st.write(f"You selected {filtered_df.shape[0]} datasets our of the {restricted_meta.shape[0]} available!")

# Convert to list and filter absolute path based on selected studies
selected_runs = [run + ".bw" for run in filtered_df["Run"].to_list()]
filtered_bigwigs = [item for item in bigwigs if any(substring in item for substring in selected_runs)]

if st.button("Next"):
    # Assign the filtered absolute paths to the analysis item
    analysis.filters.selected_studies = filtered_bigwigs
    analysis.filters.filtering_criteria = filters

    if analysis.genomic_mode:
        st.switch_page("pages/page_1.py")

    # add option to dynamically switch to page2 instead if tx mode is on
    # if analysis.TX
