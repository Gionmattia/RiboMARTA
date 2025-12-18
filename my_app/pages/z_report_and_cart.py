# z_report_and_cart.py
import streamlit as st
from module2 import *
import pandas as pd
from TranscriptID_results_helpers import *
from app_helpers import srr_lookup_sidebar

# Then I need a function to report the info as a table/something that can be exported
st.write("# Analysis Report")
st.write("Please select which analysis you want the report for")
# Choose the analysis run to visualise
chosen_run = st.pills(
    "Select the analsys run you want the report for",
    options=st.session_state.analysis_runs.keys(),
    selection_mode="single",
)
st.write(
    "Your selected option: "
    f"{ chosen_run }"
)

run = st.session_state.analysis_runs[chosen_run]

#st.stop()


# Compile the report
report_title = "RiboMARTA Report"

# Extract metadata filters and compiles neatly
filters_report = {f"{Factor}: {Levels}" for Factor, Levels in run.filters.filtering_criteria.items()}
if filters_report.__len__() > 0:
    filters_report = str(filters_report).replace("{", "").replace("}", "").replace("[","").replace("]","").replace('"', '')
elif filters_report.__len__() == 0:
    filters_report = "No filters adopted"

srr_list = {srr.split("/")[-1].replace(".bw", "") for srr in run.filters.selected_studies}
selected_studies = ", ".join(srr_list)


# Check whether analysis was genomic or transcriptomic and stores info
if run.genomic_mode.active == True:
    mode = "Genomic Coordinates"
    inputs_recap = f"""
    Mode: {mode}
    
    Coordinates selected:
        baseline: {run.genomic_mode.genomic_ranges[0]}
        target: {run.genomic_mode.genomic_ranges[1]}
        noise region: {run.genomic_mode.genomic_ranges[2]}
    
    Parameters:
        ci: {run.genomic_mode.ci}
        N: {run.genomic_mode.N}
        n_jobs: {run.genomic_mode.n_jobs}
    """

elif run.transcript_mode.active == True:
    mode = "Ensembl Transcript ID"
    inputs_recap = f"""
    Mode: {mode}
    
    Transcript ID: {run.transcript_mode.txid}

    Transcript regions selected:
        baseline: {run.transcript_mode.x_slice}
        target: {run.transcript_mode.y_slice}
        noise region: {run.transcript_mode.n_slice}

    Parameters:
        ci: {run.transcript_mode.ci}
        N: {run.transcript_mode.N}
        n_jobs: {run.transcript_mode.n_jobs}
    """

# Constructing the markdown string
markdown_content = report_content = f"""
## {report_title}
#### Analysis ID: 
{run.quirky_id}
#### Organism Selected: 
{run.organism}
#### Metadata Filters adopted: 
{filters_report}
\n Selected studies: {len(run.filters.selected_studies)} 
#### Inputs Recap:
{inputs_recap}
#### Datasets Processed (SRR):
{selected_studies}
"""

st.write(markdown_content)
print("Markdown report saved as report.md")