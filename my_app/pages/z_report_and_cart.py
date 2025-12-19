# z_report_and_cart.py
import streamlit as st
from module2 import *
import pandas as pd
from TranscriptID_results_helpers import *
from app_helpers import srr_lookup_sidebar
import io
import zipfile
from bokeh.embed import file_html
from bokeh.resources import CDN

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

if chosen_run is None:
    st.stop()

run = st.session_state.analysis_runs[chosen_run]


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

    results = run.genomic_mode.results_meta

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

    results = run.transcript_mode.txid_results_meta

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
#### Datasets Flagged:
{run.flagged_datasets}
"""

st.write(markdown_content)

# All plots are stored within the run.user_plots list (overall variable).
# The report is written within the markdown_content variable (assembled here).
# the table with the results is stored within run.transcript_mode.results_meta (mode-specific, for genomic is genomic_mode.results_meta)

def create_zip_report(markdown_content, user_plots, results):

# Create a buffer to hold the zip file
    buffy = io.BytesIO()
    with zipfile.ZipFile(buffy, "w") as zf:
        # 1) save the markdown report
        zf.writestr("report.md", markdown_content)
        # 2) save all plots as html
        for i, fig in enumerate(user_plots):
            # Convert Plotly figure to HTML string
            # Check if it's a Plotly figure
            if hasattr(fig, "to_html"):
                plot_html = fig.to_html(full_html=False, include_plotlyjs='cdn')

            # Check if it's a Bokeh figure
            else:
                try:
                    plot_html = file_html(fig, CDN, "Transcript_plot_with_ORF")
                except Exception as e:
                    plot_html = f"<h1>Error saving plot</h1><p>{str(e)}</p>"
            zf.writestr(f"plots/plot_{i + 1}.html", plot_html)
        # 3) Save results
        csv_results = results.to_csv(index=False)
        zf.writestr("results.csv", csv_results)

        return buffy


if st.button("Prepare Download Bundle"):
    zip_data = create_zip_report(markdown_content, run.user_plots, results)

    st.download_button(
        label="Download All Files (.zip)",
        data=zip_data,
        file_name="analysis_results.zip",
        mime="application/zip"
    )

