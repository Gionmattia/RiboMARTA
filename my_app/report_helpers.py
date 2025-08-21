"""
report_helpers.py

This module contains all the functions used to compile a report at the end of the analysis.
It parses through the options (and flags) the user set, to compile a human-readable, full report.

Author: Gionmattia Carancini
Date: 21-08-2025
"""

from app_helpers import *

def compile_report(current_run):

    # Get analysis mode
    analysis_mode_flags = [current_run.genomic_coordinates,
                           current_run.ensembl_transcript_id,
                           current_run.uploaded_file]
    flags_names = ["Genomic Coordinates", "Ensembl Transcript ID", "Uploaded File"]
    analysis_mode = [name for name, value in zip(flags_names, analysis_mode_flags) if value]


    # Wraps in text
    text_opener = f"""
    Analysis ID: {current_run.quirky_id}
    Organism: {current_run.organism}
    Analysis Mode: {analysis_mode}
    """

    if analysis_mode == "Genomic Coordinates":
        # Should I create a specific wrapper? like subfunctions?
        user_input_genco_ranges = current_run.genomic_ranges

    elif analysis_mode == "Ensembl Transcript ID":
        user_input_transcript_id = current_run.txid

    return text_opener


