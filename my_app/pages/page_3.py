# page_3.py
import streamlit as st
import pandas as pd
from module0 import *
from transcriptID_mode_helpers import *

st.title("File Upload Mode")

st.write("This is more of a 'last-second idea'. I thought it could be cool, but it needs more thought."
         "Overall... it's not that necessary right now, so I will just set this aside.")

st.write("At the present time, this page should update based on the state of the st.session_state.genomic_upload/txid_upload"
         "Depending on what has been offered, it could load something different on here..."
         ""
         "OR we could just redirect the inputs to page 4/5 where the analysis for genomic/transcriptomic is done and gg.."
         ""
         "")

# So there should be the option to choose one over two options... "Genomic" or "Transcript ID".

df = pd.DataFrame({
  'options': ["None Selected", "Genomic Coordinates", "Ensembl Transcript ID"]
})

option = st.selectbox(
    'Select the type of analysis',
     df["options"])

if option == 'Genomic Coordinates':
    st.write("Not implemented yet, sorry")
    st.stop()

elif option == 'Ensembl Transcript ID':

    # if Picture has not been created, creates it
    if 'picture' not in st.session_state:
        # Creates GA_df
        ga_df = global_aggregate(st.session_state.all_studies_df)
        # Makes plot
        st.session_state.picture = plot_global_aggregate_with_zero(ga_df,
                                                               st.session_state.utr5_length,
                                                               st.session_state.cds_length)
    # if picture already exists, renders it
    if st.session_state.get("picture") is not None:
        st.plotly_chart(st.session_state['picture'], use_container_width=True)

    if "show_inputs" not in st.session_state:
        st.session_state.show_inputs = False

    if st.button("Select Ranges"):
        st.session_state.show_inputs = True

    if st.session_state.show_inputs:
        col1, col2, col3 = st.columns([5, 5, 5])  # originally [5, 1, 5]
        full_mRNA_length = st.session_state.utr5_length + st.session_state.cds_length + st.session_state.utr3_length

        with col1:
            st.write("Select First Region")
            st.session_state.x_slice = ranges_selector(full_mRNA_length, "x")

        with col2:
            st.write("Select Second Region")
            st.session_state.y_slice = ranges_selector(full_mRNA_length, "y")

        with col3:
            st.write("Select Noise Region")
            st.session_state.n_slice = ranges_selector(full_mRNA_length, "n")

        if st.button("Next"):
            # NEED to add an internal check here
            st.session_state.coordinates = {"x_slice": st.session_state.x_slice,
                                            "y_slice": st.session_state.y_slice,
                                            "n_slice": st.session_state.n_slice}
            st.switch_page("pages/page_5.py")
