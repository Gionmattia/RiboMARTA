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
from mappability_helpers import bokeh_mappability_track
from mappability_helpers import txid_map_extraction
from bokeh.models import FactorRange
from bokeh.embed import file_html
from bokeh.resources import CDN
from bokeh.layouts import column
from streamlit_bokeh import streamlit_bokeh
import pandas as pd
import numpy as np
import warnings
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
from bokeh.models import BoxAnnotation, Label
from bokeh.plotting import figure, show
from bokeh.models import (
    ColumnDataSource,
    LinearColorMapper,
    ColorBar,
    HoverTool,
    BoxAnnotation,
    Label)



def txid_mode_intro():
    """
    Renders Marta and a brief description of txid mode.
    :return:
    """
    col1, col2 = st.columns([4, 2])
    with col1:
        st.markdown("Welcome to the Transcript ID Analysis page.")
        st.markdown("""1) By clicking on "Load Transcript", you will be shown a plot with the aggregate inferred A-sites of 
        the Ribosome-Protected Footprints (RPFs), across the available studies. You are also allowed to download the data as a .csv, if you wish so.
        
        **BEWARE** The positions on the transcript plots are shown as 0-indexed (the first nt positon is numbered 0, 
        instead of 1) for internal consistency. 
        
        """)
        st.markdown(
            """2) Once the plots are visualised, you can click on 'Select Ranges' to specify the areas of the transcript
             you want to compare in the analysis. You will be asked to select a baseline, a target, and a noise region.""")
        st.markdown(
            """3) Finally, you can decide whether you want to take into consideration the regions' mappability track (the dotted line
             above the plots) through the use of a Smoothing Option.""")
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
            if st.button("Load Transcript"): # If not in tx upload, loads
                current_run.loading_tx_cov = True
                # and also loads mappability track
                if analysis.mappability_path is not None:
                    mappability_df = txid_map_extraction(txid, analysis.mappability_path, st.session_state.gr)
                    mappability_df["X"] = mappability_df.index
                    other_col = [c for c in mappability_df.columns if c not in ['Position', 'X']][0]
                    current_run.mappability_df = mappability_df.rename(columns={other_col: 'Mappability'})
        elif analysis.tx_upload == True: # if tx upload mode, stores a backup
            current_run.all_studies_df_backup = current_run.all_studies_df.copy()


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


def ga_plots_UI_render(current_run, analysis):
    """
    Renders the plots stored (tx_ga_plot and triplet_ga_plot) as a st.expander.
    Also adds some info text.
    :param current_run: the class AnalysisStatus instance for the current analysis
    :return:
    """

    tab2, tab3, tab4, tab5 = st.tabs(["Raw Global Aggregate (smoothed)", "Raw Global Aggregate (unsmoothed)",
                                      "Noise Region Selector", "Normalised Global Aggregate"])
    with tab2:
        st.write("####  Raw Global Aggregate & ORF plot")
        st.write("Below you can see a Single Transcript plot for your data.")
        st.write("The x-axis represents the position on the transcript, while the y-axis represents the Global Aggregate"
                 "(the sum of raw counts from each study, in that position), smoothed by a window of size k for easy exploration.")
        st.write("Counts have been coloured by frame, and beneath the plot an ORF plot displaying starts and stops "
                 "has been added.")
        if current_run.triplet_ga_plot is not None:
            streamlit_bokeh(current_run.triplet_ga_plot, use_container_width=False)
            if st.button("Store in the report", key="triplet_ga_plot"):
                analysis.user_plots.append(current_run.triplet_ga_plot)
                st.success(f"Plot added! You currently have {len(analysis.user_plots)} plots.")

    with tab3:
        st.write("####  Raw Global Aggregate & ORF plot")
        st.write("Below you can see a Single Transcript plot for your data.")
        st.write("The x-axis represents the position on the transcript, while the y-axis represents the Global Aggregate"
                 "(the sum of raw counts from each study, in that position). Values are kept raw and without smoothing,"
                 "to allow easier visualisation of periodicity patterns when zooming in.")
        st.write("Counts have been coloured by frame, and beneath the plot an ORF plot displaying starts and stops "
                 "has been added.")
        if current_run.triplet_ga_plot_raw is not None:
            streamlit_bokeh(current_run.triplet_ga_plot_raw, use_container_width=False)
            if st.button("Store in the report", key="triplet_ga_plot_raw"):
                analysis.user_plots.append(current_run.triplet_ga_plot_raw)
                st.success(f"Plot added! You currently have {len(analysis.user_plots)} plots.")

    with tab4:
        st.write("#### Noise Control Region Selector")
        st.write("In this panel you can visualise the 3' trailer (if present) of your transcript."
        " For each position on the transcript, you can see how many of the high-counts studies share at least a mapped footprint there."
        "If a position is  conserved (= shared across studies) it will have a higher percentage."
        "\n"
        "\n Given that Noise is random, a good Noise Control Region is a section of low percentage (where mapped reads are shared by few studies)."
        "\n"
        "\n Use the green section in the plot as a reference."
                 )
        with st.expander("Options and help"):

            st.markdown(""" 
            What we are doing is: we take this datasets, and we check how many of them share at least one mapped read in their 3'UTRs.


            **Noise is random:** we expect the least conserved positions (low percentage) to be good candidates for the Noise Control region.
            You are definitely sure there is no downstream translation/pausing/anything there!
            
            ##### **IMPORTANT:** 
            This method relies on the use of only the higher-coverage datasets (too sparse data might 
            alter the plot and lead to inaccurate estimate of noise). It is important to select enough datasets to be representative, 
            but also to not include datasets with too low counts!
            
            By default, the plot requires that at least 20% of all available datasets to have high total counts (good compromise).
            
            If you get a warning, you can either:
            
            - a) Try reduce the percentage from 20% to less on the slider below (the method will work with less datasets, but still high counts).
            - b) Ignore the warning, but be aware that this might make regions look LESS conserved than what they truly are, due to low coverage.
            
            """)
                

            cutoff_text = "Select the minimum datasets to be considered for the plot. I.e. the percentage of high-coverage datasets out of the total"

            cutoff = st.slider(cutoff_text,
                          1, 60, value = 20, key =  "a_slider")

        if current_run.utr3_length is not None or current_run.utr3_length != 0:



            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                percentage_sharing, studies = density_calculator(dataframe = current_run.all_studies_df,
                                                          utr5 = current_run.utr5_length,
                                                          cds = current_run.cds_length,
                                                          cutoff = cutoff)

                for w in caught:
                    st.warning(str(w.message))

            st.write("")
            st.write(f" You are using the {studies} higher-coverage datasets out of {current_run.all_studies_df.shape[1]-1} for this plot.")
            trailer_plot = trailer_signal_plot(percentage_sharing, studies)
            streamlit_bokeh(trailer_plot, use_container_width=True, theme="light_minimal")

    with tab5:
        st.write("#### Normalised Global Aggregate & ORF plot")
        st.write("The x-axis represents the position on the transcript, while the y-axis represents the Normalised Global Aggregate"
                 "(each count in a position was normalised by the total counts of its study. The the normalised values have been summed to generate"
                 "the Normalised Global Aggregate), smoothed by a window of size k for easy exploration."
                 "The normalisation treats all studies as equivalent, regardless of the number of total counts."
                 "(ie. low-count studies contribute to the plot in the same way as high-count studies).")
        st.write("Counts have been coloured by frame, and beneath the plot an ORF plot displaying starts and stops "
                 "has been added.")
        if current_run.norm_triplet_ga_plot is not None:
            streamlit_bokeh(current_run.norm_triplet_ga_plot, use_container_width=False)


def density_calculator(dataframe, utr5, cds, cutoff = 20):

    dataset_cutoff = 100/cutoff
    # formatting
    input_df = dataframe.copy()
    if "Sequence" in input_df.columns:
        input_df = input_df.drop(columns = ["Sequence"])
    if "Position" in input_df.columns:
        input_df = input_df.drop(columns = ["Position"])

    # extract cds and trailer
    cds_df = input_df.iloc[utr5:utr5+cds]
    trailer_temp= input_df.iloc[utr5+cds:]

    # Filter studies with too low coverage (low tot counts)
    tot_counts_cds = cds_df.sum(axis = 0)
    mean = tot_counts_cds.mean()
    median = tot_counts_cds.median()
    threshold = max(median, mean)
    dataset_cutoff = len(tot_counts_cds)/dataset_cutoff # At least 1/5 of the total datasets must support this (default).

    n = (tot_counts_cds > threshold).sum()
    if n > dataset_cutoff:
        passed_indexes = tot_counts_cds[tot_counts_cds > threshold].index.tolist()
    else:
        warnings.warn("There wess less datasets among the high-counts than the percentage selected. Using all datasets for the plot") # SHOULD RAISE A WARNING
        passed_indexes = tot_counts_cds[tot_counts_cds > threshold].index.tolist()

    # filter trailer df - remove studies with too low total counts
    trailer = trailer_temp[passed_indexes].copy()

    # Encode trailer as 0 (no footprint) and 1 (footprint)
    trailer_encoded = (trailer !=0).astype(int)

    density = trailer_encoded.mean(axis=1)
    percentage_sharing = density * 100

    # report number studies used for the mean
    studies = trailer_encoded.shape[1]

    return percentage_sharing, studies


def trailer_signal_plot(series, studies):
    """
    Bokeh plot of percentage of studies sharing at least one RFP over position.
    If a position is less shared (low percentage) is more likley to be noise, rather than translation.
    """

    df_plot = series.reset_index()
    df_plot.columns = ["position", "percentage"]
    df_plot = df_plot.dropna(subset=["percentage"])

    suggested_low = df_plot["percentage"].where(df_plot["percentage"] > 0).min() * 2

    src = ColumnDataSource(df_plot)
    mapper = LinearColorMapper(palette="Viridis256", low=0.0, high=100.0)

    p = figure(
        title=f"Per-position study-contribution — low = present in few studies out of {studies} -  more likely to be noise",
        x_axis_label="3' Trailer position",
        y_axis_label="Percentage of considered datasets sharing it",
        height=380, width=950,
        tools="pan,box_zoom,wheel_zoom,reset,save",
        y_range=(0, 100),
    )

    # invisible full-height hit area — one bar per position, always hoverable
    hit = p.vbar(x="position", top=100, source=src, width=0.9,
                 fill_alpha=0, line_alpha=0)

    r = p.vbar(x="position", top="percentage", source=src, width=0.9,
               fill_color={"field": "percentage", "transform": mapper},
               line_color=None)

    p.add_layout(ColorBar(color_mapper=mapper, title="percentage"), "right")
    p.add_tools(HoverTool(renderers=[hit],
                          tooltips=[("position", "@position"),
                                    ("percentage", "@percentage{0.000}")]))

    # advisory band, drawn but clearly non-binding

    p.add_layout(BoxAnnotation(top=suggested_low, fill_alpha=0.40, fill_color="#57d677"))
    p.add_layout(Label(x=10, y=10, x_units="screen", y_units="screen",
                       text_font_size="9pt", text_color="#a0a0a0"))

    return p