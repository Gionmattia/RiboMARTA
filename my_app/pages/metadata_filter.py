import streamlit as st
import pandas as pd
from pandas.api.types import is_numeric_dtype
import os

analysis = st.session_state.analysis_runs[st.session_state.quirky_id]
st.write(f"######  :gray[ Analysis ID: {analysis.quirky_id}]")
st.title("Dataset Selection")

if analysis.genomic_upload == True:
    st.write("You are currently working with loaded data. Metadata might no be available for all your studies.")
    analysis.full_dataset_list = analysis.genomic_mode.genomic_merged.columns.tolist()
    analysis.full_dataset_list.remove("Sequence")
    analysis.full_dataset_list.remove("Position")
    restricted_meta = st.session_state.metadata[st.session_state.metadata["Run"].isin(analysis.full_dataset_list)]
    available_meta = restricted_meta.columns
    selected_columns = st.multiselect(f"Select the metadata fields to filter by", available_meta)

elif analysis.tx_upload == True:
    analysis.full_dataset_list = analysis.transcript_mode.all_studies_df.columns.tolist()
    analysis.full_dataset_list.remove("Sequence")
    analysis.full_dataset_list.remove("Position")
    st.write("You are currently working with loaded data. Metadata might no be available for all your studies.")
    restricted_meta = st.session_state.metadata[st.session_state.metadata["Run"].isin(analysis.full_dataset_list)]
    available_meta = restricted_meta.columns
    selected_columns = st.multiselect(f"Select the metadata fields to filter by", available_meta)

elif analysis.genomic_upload == False and analysis.tx_upload == False:
    # Move the dataset list to the current analysis run
    analysis.full_dataset_list = st.session_state.full_dataset_list
    # Retrieve absolute paths for all bigwigs
    bigwigs = st.session_state.bigwig_files
    # Subset the metadata based on the availability of bigwigs.
    restricted_meta = st.session_state.metadata[st.session_state.metadata["Run"].isin(analysis.full_dataset_list)]
    # Create selector based on available metadata
    available_meta =  restricted_meta.columns
    selected_columns = st.multiselect(f"Select the metadata fields to filter by", available_meta)

st.write("Here you can select which datasets to process, based on the available metadata. "
             "The current interface allows to choose from any metadata information available on the RDP,"
             " good suggestions to start from are TISSUE, CONDITION, or CELL_LINE")

# Build filters dictionary dynamically

filters = {}
for col in selected_columns:
    data = restricted_meta[col]

    # CASE 1: Continuous Data (Integers and Floats)
    if is_numeric_dtype(data) and data.nunique() > 10:
        min_val = float(data.min())
        max_val = float(data.max())

        selected_range = st.slider(
            f"Select range for {col}",
            min_value=min_val,
            max_value=max_val,
            value=(min_val, max_val)
        )
        # Only add to filters if the range has been narrowed
        if selected_range != (min_val, max_val):
            filters[col] = selected_range

    # CASE 2: Categorical Data (Strings, Booleans)
    else:
        unique_vals = data.unique()
        selected = st.multiselect(f"Select values for {col}", unique_vals)
        if selected:
            filters[col] = selected

# Apply all filtering criteria
mask = pd.Series(True, index=restricted_meta.index)

for col, values in filters.items():
    # This works if one feature is a range
    if isinstance(values, tuple):
        mask &= restricted_meta[col].between(values[0], values[1]) # &= operator acts as a "narrowing" filter
    else:
        mask &= restricted_meta[col].isin(values)

filtered_df = restricted_meta[mask]

# Message to user
st.write(f"You selected {filtered_df.shape[0]} datasets our of the {restricted_meta.shape[0]} available!")

# Convert to list and filter absolute path based on selected studies
if analysis.genomic_upload == False and analysis.tx_upload == False:
    selected_runs = [run + ".bw" for run in filtered_df["Run"].to_list()]
    filtered_bigwigs = [item for item in bigwigs if any(substring in item for substring in selected_runs)]
elif analysis.genomic_upload == True or analysis.tx_upload == True:
    selected_runs = [run + ".bw" for run in filtered_df["Run"].to_list()]
    filtered_bigwigs = selected_runs

if st.button("Next"):
    # Assign the filtered absolute paths to the analysis item
    analysis.filters.selected_studies = filtered_bigwigs
    analysis.filters.filtering_criteria = filters

    if analysis.genomic_mode.active == True:
        st.switch_page("pages/page_1.py")

    elif analysis.transcript_mode.active == True:
        st.switch_page("pages/page_2.py")