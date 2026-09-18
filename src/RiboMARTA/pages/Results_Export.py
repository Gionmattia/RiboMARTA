# z_report_and_cart.py
import streamlit as st
from module2 import *
import pandas as pd
from TranscriptID_results_helpers import *
from app_helpers import quirky_id_checkpoint
import io
import zipfile
from bokeh.embed import file_html
from bokeh.resources import CDN

quirky_id_checkpoint()
# Then I need a function to report the info as a table/something that can be exported
st.write("# Analysis Report")
st.write("Please select which analysis you want the report for")
# Choose the analysis run to visualise
chosen_run = st.pills(
    "Select the analysis run you want the report for:",
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
if run.filters.filtering_criteria is not None:
    filters_report = {f"{Factor}: {Levels}" for Factor, Levels in run.filters.filtering_criteria.items()}
    if filters_report.__len__() > 0:
        filters_report = str(filters_report).replace("{", "").replace("}", "").replace("[","").replace("]","").replace('"', '')
    elif filters_report.__len__() == 0:
        filters_report = "No filters adopted"
elif run.filters.filtering_criteria is None:
    filters_report = "No filters adopted"

if run.filters.selected_studies is not None:
    srr_list = {srr.split("/")[-1].replace(".bw", "") for srr in run.filters.selected_studies}
    selected_studies = ", ".join(srr_list)

# Check whether analysis was genomic or transcriptomic and stores info
if run.genomic_mode.active == True:
    mode = "Genomic Coordinates"

    if run.genomic_mode.smoothing_option == "Criminal Smoothing (Experimental)":
        z_score_threshold = run.genomic_mode.smoothing_criminal_thresh
    else:
        z_score_threshold =  "No z-score approach used to smooth the data (Smooth Criminal)"


    inputs_recap = f"""
    Mode: {mode}
    
    Coordinates selected:
        baseline: {run.genomic_mode.genomic_ranges[0]}
        target: {run.genomic_mode.genomic_ranges[1]}
        noise region: {run.genomic_mode.genomic_ranges[2]}
    
    Parameters:
        alpha_baseline: {run.genomic_mode.alpha_baseline}
        alpha_target: {run.genomic_mode.alpha_target}
        ci: {run.genomic_mode.ci}
        N: {run.genomic_mode.N}
        n_jobs: {run.genomic_mode.n_jobs}
        seed: {run.genomic_mode.seed}
        
    Smoothing:
        Smoothing option: {run.genomic_mode.smoothing_option}
        z-score threshold (Smooth Criminal only): {z_score_threshold}
    """

    results = run.genomic_mode.results_meta

elif run.transcript_mode.active == True:
    mode = "Ensembl Transcript ID"
    if run.transcript_mode.smoothing_option == "Criminal Smoothing (Experimental)":
        z_score_threshold = run.transcript_mode.smoothing_criminal_thresh
    else:
        z_score_threshold =  "No z-score approach used to smooth the data (Smooth Criminal)"

    inputs_recap = f"""
    Mode: {mode}
    
    Transcript ID: {run.transcript_mode.txid}

    Transcript regions selected:
        baseline: {run.transcript_mode.x_slice}
        target: {run.transcript_mode.y_slice}
        noise region: {run.transcript_mode.n_slice}
        
    Smoothing:
        Smoothing option: {run.transcript_mode.smoothing_option}
        z-score threshold (Smooth Criminal only): {z_score_threshold}
        
    Dataset Aggregation:
        Aggregation option: {run.transcript_mode.aggregation_option}
        Metadata field used for aggregation: {run.transcript_mode.meta_agg_option}

    Parameters:
        alpha_baseline: {run.transcript_mode.alpha_baseline}
        alpha_target: {run.transcript_mode.alpha_target}
        ci: {run.transcript_mode.ci}
        N: {run.transcript_mode.N}
        n_jobs: {run.transcript_mode.n_jobs}
        seed: {run.transcript_mode.seed}
    """

    results = run.transcript_mode.txid_results_meta


# Constructing the markdown string
markdown_content = report_content = f"""
## {report_title}
#### Analysis ID: 
{run.quirky_id}
#### Organism Selected: 
{run.organism}
#### Mappability track available:
{run.transcript_mode.mappability_df is not None}
#### Metadata Filters adopted: 
{filters_report}
\n Selected studies: {len(srr_list)} 
#### Inputs Recap:
{inputs_recap}
#### Datasets Processed (SRR):
{selected_studies}
#### Datasets Flagged:
{run.flagged_datasets if len(run.flagged_datasets) != 0 else "No Datasets flagged"}
"""

st.write(markdown_content)

# All plots are stored within the run.user_plots list (overall variable).
# The report is written within the markdown_content variable (assembled here).
# the table with the results is stored within run.transcript_mode.results_meta (mode-specific, for genomic is genomic_mode.results_meta)


# Temporary - I need to understand how much this works first.. but having it wrapped makes for an easier interpretation
import copy
from bokeh.models import Plot


def create_zip_report(markdown_content, user_plots, raw_data, results):
# Create a buffer to hold the zip file
    buffy = io.BytesIO()

    with zipfile.ZipFile(buffy, "w") as zf:
        # 1) save the markdown report
        zf.writestr("report.md", markdown_content)
        # 2) save all plots as html but .svg downloadable
        for i, fig in enumerate(user_plots):
            # Check if it's a Plotly figure
            # Convert Plotly figure to HTML string
            if hasattr(fig, "to_html"):
                plot_html = fig.to_html(full_html=False, include_plotlyjs='cdn')
                try:
                    zf.writestr(f"plots/plot_{i + 1}.svg", fig.to_image(format="svg"))
                except Exception as e:
                    zf.writestr(f"plots/plot_{i + 1}_svg_error.txt", str(e))
            # Check if it's a Bokeh figure
            else:
                try:
                    # We are changing the backend to svg, then back. This should avoid issues arising with
                    # copy() method (not working on plots) or shallow copying that could break shaered children.
                    plots = list(fig.select({'type': Plot}))
                    previous = [p.output_backend for p in plots]
                    for p in plots:
                        p.output_backend = "svg"
                    html_svg = file_html(fig, CDN, f"Plot {i + 1}")
                    for p, prev in zip(plots, previous):
                        p.output_backend = prev
                    zf.writestr(f"plots/plot_{i + 1}.html", html_svg)
                except Exception as e:
                    zf.writestr(f"plots/plot_{i + 1}_svg_error.txt", str(e))

            #zf.writestr(f"plots/plot_{i + 1}.html", plot_html)
        # 3) Store raw data (for reference) NB. does not store normalised/aggregated data! just RAW as it was extracted
        csv_raw_data = raw_data.to_csv(index=True)
        zf.writestr("raw_data.csv", csv_raw_data)
        # 4) Save results
        if results is not None:
            csv_results = results.to_csv(index=True)
            zf.writestr("results.csv", csv_results)

        return buffy

col1, col2, col3 = st.columns([2, 2, 2])
with col1:
    filename = st.text_input("Select name for download file", "analysis_results")
    filename = f"{filename}.zip"

    st.write(f"Your results will be saved as {filename}")

if st.button("Prepare Download Bundle"):
    zip_data = create_zip_report(markdown_content = markdown_content,
                                 user_plots = run.user_plots,
                                 raw_data= run.transcript_mode.all_studies_df_backup,
                                 results =  results)

    st.download_button(
        label="Download All Files (.zip)",
        data=zip_data,
        file_name=filename,
        mime="application/zip"
    )

