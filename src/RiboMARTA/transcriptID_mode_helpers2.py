"""
transcriptID_mode_helpers2.py

This module contains functions for the transcriptID mode. Specifically, it contains functions used to render text, plots,
and other "front-end" options.
For functions relative to backend operations, such as extracting transcript_ranges or generating the plots in the first place,
please refer to transcriptID_mode_helpers.py

Author: Gionmattia Carancini
Date: 25-08-2025
"""
import streamlit as st
from app_helpers import AnalysisStatus
from transcriptID_mode_helpers import *
from mappability_helpers import bokeh_mappability_track # Needed
from mappability_helpers import txid_map_extraction # Needed
# from bokeh.models import FactorRange
# from bokeh.embed import file_html
# from bokeh.resources import CDN
# from bokeh.layouts import column
from streamlit_bokeh import streamlit_bokeh
import pandas as pd
import numpy as np
import warnings
# from statsmodels.stats.multitest import multipletests
# from joblib import Parallel, delayed
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
    Renders Marta and a brief description of txid mode
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
        st.image(st.session_state["marta_images"]["explaining"], width=180)


def txid_checkpoint(current_run):
    """
    Checkpoint function. It checks if the user has submitted a transcriptID before allowing
    to further proceed in the analysis.
    Txid should be saved in the current_run (class: AnalysisStatus, see app_helpers.py)
    """

    txid = current_run.txid
    if txid == "":
        col1, col2 = st.columns([1, 3])
        with col1:
            st.image(st.session_state["marta_images"]["error"], width=127)
        with col2:
            st.markdown("<br><br>", unsafe_allow_html=True)
            st.error(" The 'Transcript ID' field was empty. Please provide a value before proceeding!")
        st.stop()

    return

def format_txid_df_download(current_run):
    """
    Sets up the df with coverage from all Runs for download.
    It adds the sequence information to it as well.
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
        file_name=f"{current_run.txid}_5'_{current_run.utr5_length}_CDS_{current_run.cds_length}_3'_{current_run.utr3_length}_raw_counts.csv",
        mime='text/csv'
    )


def ga_plots_UI_render(current_run, analysis):
    """
    Renders the plots stored (tx_ga_plot and triplet_ga_plot) as a st.expander.
    Also adds some info text.
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
        st.write("In this panel you can visualise the 3' trailer of the transcript (or the whole, if not CDS is annotated)."
        " For each position on the transcript, you can see how many of the high-coverage studies share at least a mapped footprint there."
        "If a position is  conserved (= shared across studies) it will have a higher percentage."
        "\n"
        "\n Given that Noise is random, a good Noise Control Region is a section of low percentage (where mapped reads are shared by few studies)."
        "\n"
        "\n Use the green section in the plot as a reference."
                 )
        with st.expander("Options and help for the scared"):

            st.markdown(""" 
            
            Ok so:

            **Noise is random:** we expect the least conserved positions (low percentage) to be good candidates for the 
            Noise Control region. You can assume there is no downstream translation/pausing/anything there!
            
            This line of reasoning holds so long as we assume to have high-coverage datasets.
            
            We first select the higher-coverage datasets. If the transcript has a CDS, we check how many total counts they have there.
            If they don't have a CDS, we use the full transcript for the plot.
            
            We consider only datasets that have tot counts above the mean for all datasets.

            Since this might lead to a few high-coverage dominating the profile, there's the cutoff slider. You decide
            what minimum percentage of the total is required.
            
            If after filtering you have less than cutoff, the plot issues a warning and defaults to use all datasets.
            You can:
            
            - a) Try to reduce the percentage from 20% to less on the slider below (the method will work with fewer datasets, but still high counts).
            - b) Ignore the warning, but be aware that this might make regions look LESS conserved than what they truly are, due to low coverage.
            
            Nb. This is a SUGGESTION, nothing beats your eye!
            """)
                

            cutoff_text = "Select the minimum datasets to be considered for the plot. I.e. the percentage of high-coverage datasets out of the total"

            cutoff = st.slider(cutoff_text,
                          1, 60, value = 20, key =  "a_slider")

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                # Failsafe - if annotation is missing for either, it just defaults to 0
                # density calculator will use the whole transcript for the OC plot
                utr5_length = current_run.utr5_length if current_run.utr5_length is not None else 0
                cds_length = current_run.cds_length if current_run.cds_length is not None else 0

                percentage_sharing, studies, n_total = density_calculator(dataframe = current_run.all_studies_df,
                                                                          utr5 = utr5_length,
                                                                          cds = cds_length,
                                                                          cutoff = cutoff)

                for w in caught:
                    st.warning(str(w.message))

            if percentage_sharing.empty:
                st.warning("The transcript isoform ended with the annotated CDS. We advise to use another one or to "
                           "resort to RiboSeq.Org GWIPS + Genomic Coordinates mode to sample noise.")
            else:
                st.write("")
                if studies == n_total:
                    st.write(f"You are using all datasets ({n_total}) for this plot")
                else:
                    st.write(f"You are using the {studies} higher-coverage datasets out of {n_total} for this plot.")

                OC_plot = trailer_signal_plot(percentage_sharing, studies)
                streamlit_bokeh(OC_plot, use_container_width=True, theme="light_minimal")

            # col1, col2 = st.columns([2, 6])
            # with col1:
            #     st.write("")
            #     # window_size = st.number_input("Select the size of the noise region you want to get", value = 100, key="n_window")
            # with col2:
            #     st.write("")
            #     st.write("")
            #     st.write("Suggested Intervals not yet implemented. Requires extra work")
            #     # give_intervals(percentage_sharing, int(window_size))

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


def density_calculator(dataframe: pd.DataFrame,
                       utr5: int,
                       cds: int,
                       cutoff: int = 20) -> tuple[pd.Series, int, int]:
    """
    Function to compute the data for the OC plot: what % of the high-coverage datasets has (at least) a coutn mapped to
     the position of the trailer?
        If no CDS annotated, then it uses whole transcript instead.
    Parameters
    -------
        dataframe: pd.DataFrame with all coverage per position for the selected transcript
        utr5: int length of the 5' leader
        cds: int length of the CDS
        cutoff: percentage of studies we want to consider as high-coverage datasets

    Returns
    -------
        tuple: series with the % of high-cov datasets with at least 1 count per position, and number of datasets
    """
    # formatting
    input_df = dataframe.copy()
    if "Sequence" in input_df.columns:
        input_df = input_df.drop(columns = ["Sequence"])
    if "Position" in input_df.columns:
        input_df = input_df.drop(columns = ["Position"])

    # We have an annotated CDS on the transcript:
    if cds > 0:
        counts_df = input_df.iloc[utr5:utr5+cds]
        trail_df_temp= input_df.iloc[utr5+cds:]
        if trail_df_temp.empty:
            warnings.warn("The CDS ends at the end of the annotated transcript: there is no 3' trailer to use for "
                          "sampling noise. We suggest using Genomic Coordinates mode instead.")
            return pd.Series(dtype=float), 0, input_df.shape[1]
    else:
        counts_df = input_df
        trail_df_temp = input_df

    tot_counts_cds = counts_df.sum(axis = 0)
    n_total = len(tot_counts_cds)
    # Select datasets where coverage is above mean
    high_coverage = tot_counts_cds[tot_counts_cds > tot_counts_cds.mean()]
    percentage_passed = len(high_coverage)/ n_total

    if percentage_passed > cutoff /100:
        passed_indexes = high_coverage.index.tolist()
        # failsafe = False
    else:
        warnings.warn(f"Only {percentage_passed:.0%} datasets have coverage above the mean, "
                      f"but the current cutoff requires at least > {cutoff}%. "
                      f"Using all datasets for the profile.")

        passed_indexes = tot_counts_cds.index.tolist()
        # failsafe = True

    # filter trail df - remove studies with too low coverage (ndr. could be 3' trailer or "everything". Trail is just a variable name)
    trail = trail_df_temp[passed_indexes].copy()

    # Encode trailer as 0 (no footprint) and 1 (footprint)
    trail_encoded = (trail !=0).astype(int)

    density = trail_encoded.mean(axis=1)

    # Output 1 - pd.Series that reports percentage of datasets per position
    percentage_sharing = density * 100

    # Output 2 - number studies used
    studies = trail_encoded.shape[1]

    return percentage_sharing, studies, n_total


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
        y_axis_label="Percentage of considered datasets sharing at least one count",
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

# Not sure I like this code. Too many silent issues,,, it was supposed to be easier, honestly.
# def give_intervals(percentage_sharing, window):
#     # Not sure of behaviour on edge-cases
#     rolling_mean = percentage_sharing.rolling(window=window, center=True).mean()
#
#     if rolling_mean.notna().sum() == 0:
#         st.write(f"Trailer too short for a {window}-nt window.")
#         return
#
#     min_pos = int(np.nanargmin(rolling_mean.values))
#
#     n = len(percentage_sharing)
#     half = window // 2
#
#     start = max(0, min_pos - half)
#     end   = min(n, start + window)
#     start = max(0, end - window)
#
#     start_label = percentage_sharing.index[start]
#     end_label = percentage_sharing.index[end - 1]
#
#     st.write(f"Suggested Interval (please check against plot): {start_label} - {end_label}") # return embedded within function. Nope
#     return # start, end