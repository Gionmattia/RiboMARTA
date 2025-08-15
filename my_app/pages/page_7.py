import streamlit as st
from chisquared_helpers import *
from TranscriptID_results_helpers import *
import pandas as pd
from breakdown_helpers import *

st.title("Chi Squared Tests (v.0.3)")

current_run = st.session_state.runs[st.session_state.quirky_id]

# Function to calculate the two metrics in a single dataframe.
st.write("## NEW IDEA 2")
st.write("In this case, we are actually taking into accounts Total counts, and a score which represents the percentage."
         "We unit-step transform our data, so it's either 0 (no footprint) or 1 (footprint)."
         "This allows us to compute a percentage for each region, that tells us 'what percentage of positions have at last one count?'"
         "now, with this breakdown, we can stratify our data. Datasets where all three regions are less than 0.05 are to be considered poor."
         "or, generally speaking, unreliable for the breakdown calculations."
         "(this should apply generally well, unless transcript have a very small ORF compared to the overall length)"
         " Datasets where we are between 0.05 and 0.1 are so and so... above they start to be good."
         "overall, we can argue that datasets with High total counts, and with a high percentage in at least one region are the best ones."
         "If the tot counts is high, but the percentage is poor, then they are not really reliable.")

raw_counts_breakdowns, totals_and_percents_breakdowns = get_breakdowns_slices(current_run.all_studies_df,
                           current_run.txid_ranges, # could also use current_run.genomic_merged_ranges
                           current_run.txid_results)

raw_counts_breakdowns
totals_and_percents_breakdowns

a, b =get_breakdowns_full_tx(current_run.all_studies_df,
                       current_run.utr5_length,
                       current_run.cds_length,
                       current_run.txid_results)

a
b
st.stop()
# Does the breakdowns (leader, cds, trailer) for every run and the Global Aggregate
df3 = get_breakdowns(current_run.all_studies_df,
                     current_run.txid_results,
                     current_run.utr5_length,
                     current_run.cds_length)


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

