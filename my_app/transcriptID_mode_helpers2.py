"""
transcriptID_mode_helpers2.py

This module contains functions for the transcriptID mode. Specifically, it contains functions used to render text, plots,
and other "front-end" options.
For functions relative to backend operations, such as extracting ranges or generating the plots in the first place,
please refer to transcriptID_mode_helpers.py

Author: Gionmattia Carancini
Date: 25-08-2025
"""
import streamlit as st
from app_helpers import AnalysisStatus
from transcriptID_mode_helpers import *
from bokeh.models import FactorRange
from bokeh.embed import file_html
from bokeh.resources import CDN
from bokeh.layouts import column


def txid_mode_intro():
    """
    Renders Marta and a brief description of txid mode.
    :return:
    """
    col1, col2 = st.columns([4, 1])
    with col1:
        st.write("Welcome to the Transcript ID Analysis page.")
        st.write("""1) By clicking on "Load Transcript", you will be shown a plot with the aggregate Ribo-Seq footprints 
        across all the available studies. You are also allowed to download the data as a .csv, if you wish so.""")
        st.write(
            """2) Once the plots are visualised, you can click on 'Select Ranges' to specify the areas of the transcript
             you want to compare in the analysis. You will be asked to select a baseline, a target, and a noise region.""")
        st.write(
            """3) Finally, you can decide whether you want to take into consideration the region's reading frames 
            (frame-aware analysis) or not.""")
    with col2:
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.image("martas/explaining_marta.png", width=180)


def txid_start_checkpoint(current_run, analysis):
    """
    Checkpoint function. It checks if the user has submitted a transcriptID before allowing
    to further proceed in the analysis.
    Txid should be saved in the current_run (class: AnalysisStatus, see app_helpers.py)
    It overwrites the value of current_run.loading_tx_cov flag tp True (used for loading screen)
    :param current_run:
    :return:
    """

    txid = current_run.txid
    if txid == "":
        col1, col2 = st.columns([1, 3])
        with col1:
            st.image("martas/error_marta.png", width=127)
        with col2:
            st.markdown("<br><br>", unsafe_allow_html=True)
            st.error(" The 'Transcript ID' field was empty. Please provide a value before proceeding!")
    # Allows to extract the tx coverage, if not in Upload Mode.
    else:
        if analysis.tx_upload == False:
            if st.button("Load Transcript"):
                current_run.loading_tx_cov = True


def format_txid_df_download(current_run):
    """
    Sets up the df with coverage from all Runs for download.
    It adds the sequence information to it as well.
    :param current_run:
    :return:
    """
    # copy the current working dataframe
    export_dataframe = current_run.all_studies_df.copy()
    # Add the sequence to the copy, and allow for download.
    export_dataframe.insert(0, 'Sequence', list(current_run.transcript_sequence))
    # Allow to save the raw data as a .csv
    csv = export_dataframe.to_csv(index=True).encode('utf-8')
    return st.download_button(
        label="Download data as CSV",
        data=csv,
        file_name=f'{current_run.txid}_raw_counts.csv',
        mime='text/csv'
    )


def plot_ga_with_frames(
        current_run: AnalysisStatus,
        ga_df: pd.DataFrame,
        y_column: str,
        smoothing: bool = False
) -> column:
    """
    Create a frame-aware global aggregate plot with ORF tracks beneath.

    Generates a multi-panel Bokeh visualization showing:
    - Top panel: Global aggregate coverage colored by reading frame
    - Bottom panels: ORF tracks for each of the three reading frames

    Args:
        current_run: Analysis state containing transcript metadata and sequence
        ga_df: DataFrame with 'Global Aggregate' column containing coverage data
        y_column: Name of the column to plot on y-axis
        smoothing: If True, smooth data by frame-specific windows; if False,
                   use raw values (default: False)

    Returns:
        Bokeh column layout with transcript plot and three ORF track panels
    """
    # prepare data: add frame info and optionally smooths the data.
    processed_df, k = smoother(ga_df, smoothing)

    # Creates shared x-axis as FactorRange object (allows for multiple plots to share X)
    positions = processed_df['X'].astype(str).tolist()
    x_range = FactorRange(*positions)

    # Create main transcript plot
    transcript_plot = subplot_single_transcript_plot(
        df = processed_df,
        x_range = x_range,
        y_axis = y_column,
        utr5_length=current_run.utr5_length,
        cds_length=current_run.cds_length,
        txid=current_run.txid,
        k=k)

    transcript_plot.height = 400  # or whatever height you need
    transcript_plot.sizing_mode = 'stretch_width'

    # Code to retrieve sequence and codon positions
    orf_info = compute_codons_by_frame(current_run.transcript_sequence)

    # Create ORF plots for each frame
    orf_plots = []
    for frame in [1,2,3]:
        orf_plot = ORF_plotter(orf_info, frame, x_range)
        orf_plot.height = 75
        orf_plot.sizing_mode = 'stretch_width'
        orf_plots.append(orf_plot)

    # Combine all plots into a single one
    combined_plot = column(transcript_plot,
                           orf_plots[0],
                           orf_plots[1],
                           orf_plots[2],
                           sizing_mode='stretch_width')

    return combined_plot


def plot_norm_ga_with_frames(current_run, ga_df):
    """
    Makes the plot with triplet periodicity info, and the ORF plots beneath it..

    :param current_run: the class AnalysisStatus instance for the current analyisis
    :param ga_df: dataframe with global aggregate

    The plot generated is saved within the current_run object, as "triplet_ga_plot"
    """
    # Calculates frame and smoothes data (smoothing ny window size k)
    smoothed_ga_df, k = smoother(ga_df)
    # Defines x-axis as FactorRange object (allows for multiple plots to share X)
    positions = smoothed_ga_df['X'].astype(str).tolist()
    x_range = FactorRange(*positions)
    picture2_transcript_plot = subplot_single_transcript_plot(smoothed_ga_df,
                                                              'Global_Aggregate_smoothed',
                                                              x_range,
                                                              k,
                                                              start=current_run.utr5_length,  # NEED TO FIX THIS
                                                              stop=current_run.cds_length + current_run.utr5_length,
                                                              txid=current_run.txid)

    picture2_transcript_plot.height = 400  # or whatever height you need
    picture2_transcript_plot.sizing_mode = 'stretch_width'

    # Code to retrieve sequence and codon positions
    orf_info = compute_codons_by_frame(current_run.transcript_sequence)
    # Code to create the ORF plots and the combined plot.
    orfplot_f1 = ORF_plotter(orf_info, 1, x_range)
    orfplot_f1.height = 75
    orfplot_f1.sizing_mode = 'stretch_width'
    orfplot_f2 = ORF_plotter(orf_info, 2, x_range)
    orfplot_f2.height = 75
    orfplot_f2.sizing_mode = 'stretch_width'
    orfplot_f3 = ORF_plotter(orf_info, 3, x_range)
    orfplot_f3.height = 75
    orfplot_f3.sizing_mode = 'stretch_width'
    plot = column(picture2_transcript_plot, orfplot_f1, orfplot_f2, orfplot_f3,
                                         sizing_mode='stretch_width')
    return plot


def ga_plots_UI_render(current_run):
    """
    Renders the plots stored (tx_ga_plot and triplet_ga_plot) as a st.expander.
    Also adds some info text.
    :param current_run: the class AnalysisStatus instance for the current analyisis
    :return:
    """
    with st.expander("Global Aggregate Plots"):
        tab1, tab2, tab3 = st.tabs(["Simple Aggregate", "Raw Global Aggregate & ORF plot", "Normalised Global Aggregate & ORF plot"])
        with tab1:
            st.write("#### Simple Global Aggregate plot")
            st.write("Below you can see a Single Transcript plot for your data.")
            st.write("The x-axis represents the position on the transcript, while the y-axis represents the Global Aggregate"
                     "(the sum of raw counts from each study, in that position).")
            st.write("This plot is a simplified visualisation, with no information on ORFs or reading frames."
                     "It's an easy-on-the eyes visualisation, less informative but useful for an initial exploration.")
            if current_run.tx_ga_plot is not None:
                st.plotly_chart(current_run.tx_ga_plot, use_container_width=True)
        with tab2:
            st.write("####  Raw Global Aggregate & ORF plot")
            st.write("Below you can see a Single Transcript plot for your data.")
            st.write("The x-axis represents the position on the transcript, while the y-axis represents the Global Aggregate"
                     "(the sum of raw counts from each study, in that position), smoothed by a window of size k for easy exploration.")
            st.write("Counts have been coloured by frame, and beneath the plot an ORF plot displaying starts and stops "
                     "has been added.")
            if current_run.triplet_ga_plot is not None:
                # Patchy, maybe I can come up with something better... but does it job for now.
                html = file_html(current_run.triplet_ga_plot, CDN, "Transcript_plot_with_ORF")
                st.components.v1.html(html, height=650, width=1600)
                #html = file_html(current_run.triplet_ga_plot, CDN, "Transcript_plot_with_ORF")
                #st.components.v1.html(html, height=650, width=1600)
        with tab3:
            st.write("#### Normalised Global Aggregate & ORF plot")
            st.write("The x-axis represents the position on the transcript, while the y-axis represents the Normalised Global Aggregate"
                     "(each count in a position was normalised by the total counts of its study. The the normalised values have been summed to generate"
                     "the Normalised Global Aggregate), smoothed by a window of size k for easy exploration."
                     "The normalisation treats all studies as equivalent, regardless of the number of total counts."
                     "(ie. low-count studies contribute to the plot in the same way as high-count studies).")
            st.write("Counts have been coloured by frame, and beneath the plot an ORF plot displaying starts and stops "
                     "has been added.")
            if current_run.norm_triplet_ga_plot is not None:
                # Saving as html, since blokeh v.3+ is not directly compatible with streamlit...
                # Patchy, maybe I can come up with something better... but does it job for now.
                html = file_html(current_run.norm_triplet_ga_plot, CDN, "Normalised_transcript_plot_with_ORF")
                st.components.v1.html(html, height=650, width=1600)


def frame_aware_analysis_intro():
    with st.popover("About the Frame-Aware Analysis",
                    use_container_width=True):

        st.write("### About Frame Awareness")
        st.write("The analysis uses coverage data between regions coming from all three different reading frames."
                 "While this allows to have stronger signal, it is less precise and cannot allow for the discrimination"
                 "of overlapping or nested ORFs.")
        st.write("The Frame-Aware analysis allows to compute the RTA using signal coming from each reading frame "
                 "individually. While this is mode is more precise, it has the downside not all individual studies will"
                 " have an in-frame signal strong enough to be detected as significant.")
        st.write("")