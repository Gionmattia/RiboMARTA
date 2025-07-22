import streamlit as st
from chisquared_helpers import *
from TranscriptID_results_helpers import *

st.title("Chi Squared Tests (alpha)")

st.write("""In this page, you can observe the results of the X-squared test of homogeneity, over the runs.
This tests aims at testing the hypothesis of association between one run and the distribution of reads across 
5' leader, CDS, and 3' trailer.
Since we are using big numbers, p-values willalways be significant. You must look at the cramer's V to have an 
idea on the effect size. if it is above 0.1, it's interesting. If it's above 0.5, that's jackpot.

If that's the case, look at the heatmap to see which runs/regions contributed the most to the. The strongest the colour, 
the greater the standard deviation and hence the contribution""")


df = get_breakdowns_percent(st.session_state.all_studies_df,
                    st.session_state.txid_results,
                    st.session_state.utr5_length,
                    st.session_state.cds_length)

chi2, p_value, dof, cramers_v, std_residuals = execute_chi_squared(df)

st.write("## Results")
st.write("Chi-squared: ", chi2)
st.write("p-value: ", p_value)
st.write("dof: ", dof)
st.write("cramer's V: ", cramers_v)

st.write("## Table of standard residuals")
st.write("This table is useful to see what regions/Runs differed the most. Useful if cramer's V above 0.1. ")
std_figure = plot_chisquared_residuals(std_residuals)

st.plotly_chart(std_figure, use_container_width=True)

run_option = st.multiselect("Select one or more (max 5) different Runs to plot:",
                            st.session_state.results_meta.index)

st.write("Manually check the RUN here")
ml_plot = multi_line_plot(st.session_state.all_studies_df,
                          st.session_state.utr5_length,
                          st.session_state.cds_length,
                          run_option,
                          st.session_state.coordinates)
st.plotly_chart(ml_plot, use_container_width=True)

st.write("## Going by TISSUE")
st.write("Select one of the metadata, do the chi-squared over it only? not sure if that's worth checking, if individual"
         "SRR do not yield any meaningful insight")

st.write("I think it could be more interesting to compare the actual slices (x, y) among the different runs..."
         "It could be a different perspective on RiboMARTA and be more informative.")


coordinates = st.session_state.get("coordinates", None)

st.write(f"""You selected the following ranges for the analysis:

- First Region - From {coordinates["x_slice"][0]} to {coordinates["x_slice"][1]}
- Second Region - From {coordinates["y_slice"][0]} to {coordinates["y_slice"][1]}
""")

df2 = get_breakdowns_slices(st.session_state.all_studies_df,
                            st.session_state.txid_results,
                            coordinates)

chi2, p_value, dof, cramers_v, std_residuals = execute_chi_squared(df2)

st.write("## Results")
st.write("Chi-squared: ", chi2)
st.write("p-value: ", p_value)
st.write("dof: ", dof)
st.write("cramer's V: ", cramers_v)

st.write("## Table of standard residuals")
st.write("This table is useful to see what regions/Runs differed the most. Useful if cramer's V above 0.1. ")
std_figure_2 = plot_chisquared_residuals(std_residuals)

st.plotly_chart(std_figure_2, use_container_width=True)


st.write("""So far it's just telling me there is a difference in the distributions of raw counts, across all these Runs.
It will also show in which specific cell I am observing this difference, if I look at the heatmap.

Problem is, it is too sensitive. Too many false positives, due to the high number of counts.

It's could be a filtering. But again, it has TOO MUCH sensitivity -> type I error (reject Null Hypothesis).


I resorted to use the "percentage of translational activity" (so the sum of counts per area, in percentage) as a measure.
This seems to work a bit better. Homogenity works fine...
 
Though I wonder if I could also do a pairwise, where I compare each individual one to the normalised global aggregate, converted to percentage.
That way I could have a "average of everything" profile, to which I compare each individual one, looking for striking differences.

It's interesting, but that also ponders the question...could I make a heatmap like the one I saw before then, showing the differences
from the average in term of proportions?

Would the test be statistically significant because it's too sensitive, or will it work properly now?

- what if every case is too different from the average -> everything is signficant... that would eb a nightmare
""")


st.write("This table is useful to see what regions/Runs differed the most. Useful if cramer's V above 0.1. ")
# st.write("""RiboMARTA - Roadmap for reference.
# - IDEAS: introduce the X-squared breakdowns. [alpha done]
# - IDEAS: introudce the effieniency method (after testing).
# - Introduce a "select organism" function that loads organisms. Load one, cache, forget about it. [done]
# - Implement auto-noise selection (picks a range of low z-scores from transcript).
# - Option to run test in a different way, if 3'UTR not present or too small (was it wilcoxon?).
# """)
