import streamlit as st
from chisquared_helpers import *
from TranscriptID_results_helpers import *
import pandas as pd
from breakdown_helpers import *
from genomic_mode_helpers import multi_line_plot_genomic

current_run = st.session_state.runs[st.session_state.quirky_id]
st.write(f"######  :gray[ Analysis ID: {st.session_state.quirky_id}]")

st.title("Chi Squared Tests  and breakdowns")

st.write("##Work in progress! 🚧")
col1, col2, col3 = st.columns([3,1,2])
with col1:
    st.write("This page is currently under development. Here you can see an early preview."
             "Future updates include:"
             ""
             "- Improvements to the analysis flow control"
             
             "- Introduction of a 'frame-aware' analysis mode (for nested ORFs)"
             
             "- Codon usage analysis, between selected regions (and possible correction if there's usage bias)."
             
             "- Big Beautiful 'Download ALL' Button")
with col3:
    st.image("martas/WIP_marta.png", width = 200)



st.write("In this page you can see if there's a difference in counts proportions across the three regions of interest."
         "All the results below are already significant (so there IS a difference) between each individual run and the"
         "Global Aggregate, if we look at proportions. The cramer's v will give you an idea on the effect size. Values below "
         "0.2 are to be considered irrelevant: the effect size is neglegible"
         "Studies are represented in the heatmap in descending order of total counts, with the first row representing"
         "the proportions for the Global Aggregate. A study with a high cramer's v at the top is to be considered more"
         "reliable than one at the bottom, where counts might too sparse to back up any hypothesis.")

# NEED TO USE FLAGS BETTER, THIS IS REALLY PROVISIONARY
if current_run.genomic_coordinates == True and current_run.regions_extracted == True:

# Code to visualise breakdowns from genomic coord mode:
    raw_counts_breakdowns, totals_and_percents_breakdowns, ga = get_breakdowns_slices(current_run.genomic_merged,
                                                                                     current_run.genomic_merged_ranges,
                                                                                     current_run.genco_results)

    results_df_raw = perform_pairwise_chi(raw_counts_breakdowns, baseline_idx = "Global Aggregate")

    df = prepare_chisquared_results(results_df_raw, totals_and_percents_breakdowns, ga)

    fig4 = better_plotter(df, ["% First Region", "% Second Region", "% Noise Region"])
    st.plotly_chart(fig4, use_container_width=True, key = "genco_breakdowns")

elif current_run.ensembl_transcript_id == True and current_run.txid_analysis_done == True:
    st.write( "In this page you can see if there's a difference in counts proportions across the three regions of interest."
        "Please, select if you are interested in observing the breakdowns over the selected transcript ranges,"
        "or if you want to see across 5' leader, CDS, and 3' trailer")

    col1, col2, col3 = st.columns([1, 1, 4])
    with col1:
        if st.button("Transcript Regions"):
            current_run.tx_regions = True
            current_run.tx_ranges = False
    with col2:
        if st.button("Transcript Ranges"):
            current_run.tx_regions = False
            current_run.tx_ranges = True

    if current_run.tx_ranges == True:

        raw_counts_breakdowns, totals_and_percents_breakdowns, ga = get_breakdowns_slices(current_run.all_studies_df,
                                                                                          current_run.txid_ranges, # could also use current_run.genomic_merged_ranges
                                                                                          current_run.txid_results)

        results_df_raw = perform_pairwise_chi(raw_counts_breakdowns, baseline_idx = "Global Aggregate")

        df = prepare_chisquared_results(results_df_raw, totals_and_percents_breakdowns, ga)


        fig4 = better_plotter(df, ["% First Region", "% Second Region", "% Noise Region"])
        st.plotly_chart(fig4, use_container_width=True, key = "breakdowns_txid_ranges")

    elif current_run.tx_regions == True:

        raw_counts_breakdowns, totals_and_percents_breakdowns, ga = get_breakdowns_full_tx(current_run.all_studies_df,
                                                                                           current_run.utr5_length,
                                                                                           current_run.cds_length,
                                                                                           current_run.txid_results)


        results_df_raw = perform_pairwise_chi(raw_counts_breakdowns, baseline_idx="Global Aggregate")

        df = prepare_chisquared_results(results_df_raw, totals_and_percents_breakdowns, ga)

        fig4 = better_plotter(df, ["% First Region", "% Second Region", "% Noise Region"])
        st.plotly_chart(fig4, use_container_width=True, key="breakdowns_txid_ranges")

else:
    st.stop()

with st.expander("Single study visualiser"):
    st.write("Here you can manually check the footprints for a study that picks your interest.")
    # Manual inspection of datasets


    if current_run.genomic_coordinates == True and current_run.regions_extracted == True:
        st.write("Manually check the RUN here")
        run_option = st.multiselect("Select one or more (max 5) different Runs to plot:",  current_run.results_meta.index)
        ml_plot = multi_line_plot_genomic(current_run.genomic_merged,
                                          run_option,
                                          current_run.genomic_merged_ranges)
        st.plotly_chart(ml_plot, use_container_width=True, key = "genco_single_run_plot")

    elif current_run.ensembl_transcript_id == True and current_run.txid_analysis_done == True:

        st.write("Manually check the RUN here")
        run_option = st.multiselect("Select one or more (max 5) different Runs to plot:",
                                    current_run.txid_results_meta.index)
        ml_plot = multi_line_plot(current_run.all_studies_df,
                                  current_run.utr5_length,
                                  current_run.cds_length,
                                  run_option,
                                  current_run.txid_ranges)
        st.plotly_chart(ml_plot, use_container_width=True, key = "txid_single_run_plot")




st.stop()
# Chi squared, done pairwise and using the breakdowns obtained from raw counts.
results_df_raw = perform_pairwise_chi(df3, baseline_idx = "Global Aggregate")

# Compute metrics: total counts per dataset, and percentage of transcript with at least 1 count
metrics_df = get_metrics(current_run.all_studies_df, current_run.txid_results)

# Combine the results to the metrics tot counts and percent cover
new_df = pd.merge(results_df_raw, metrics_df, on = "Run")


percentage_breaks = get_breakdowns_percent_with_GA(
    current_run.all_studies_df,
    current_run.txid_results,
    current_run.utr5_length,
    current_run.cds_length)

perc_GA = percentage_breaks[percentage_breaks["Run"]=="Global Aggregate"]


final_df = pd.merge(percentage_breaks, new_df, on = "Run")

final_df

# plot the p-values as a result
fig3_5 = pairwise_comparison_plotter(results_df_raw)
st.plotly_chart(fig3_5, use_container_width=True, key = "cmon")


# Need to create the new heatmap function.... we are almost done.
with st.expander("Detailed Breakdown Plot"):
    st.write("Here you can see a plot displaying if there is any difference in density across the three regions of the"
             "transcript. While the plot displays all significant changes, the 'cramer v' value informs of the effect size"
             "a v less than 0.1 indicates the runis neglegible, as it means that, despite being significant, the effect size is little to none."
             "Studies with a higher total counts (and hence more reliable) are displayed top, with numbers descending")

    cleaned_df = final_df[final_df['bonferroni_p'] < 0.05].copy()
    cleaned_df = cleaned_df.sort_values(by="tot_counts", ascending=True)
    # Add back the Global Aggregate
    perc_GA = percentage_breaks[percentage_breaks["Run"] == "Global Aggregate"]
    values = df3[df3["Run"] == "Global Aggregate" ].drop(columns="Run")
    tot = values.sum(axis=1)
    GA_vals = {"UTR_5": perc_GA["UTR_5"][0],
               "UTR_3": perc_GA["UTR_3"][0],
               "CDS": perc_GA["CDS"][0],
               "Run": "Global Aggregate",
               "tot_counts": int(tot),
               "cramers_v": 0,
               "bonferroni_p": 1}
    cleaned_df.loc[len(cleaned_df)] = [GA_vals.get(c, None) for c in cleaned_df.columns]

    fig4 = better_plotter(cleaned_df)
    st.plotly_chart(fig4, use_container_width=True)


# Manual inspection of datasets
run_option = st.multiselect("Select one or more (max 5) different Runs to plot:",
                            current_run.txid_results_meta.index)

st.write("Manually check the RUN here")
ml_plot = multi_line_plot_simple(current_run.all_studies_df,
                          current_run.utr5_length,
                          current_run.cds_length,
                          run_option)
st.plotly_chart(ml_plot, use_container_width=True)

