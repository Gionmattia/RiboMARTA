import streamlit as st
from transcriptID_mode_helpers import *

st.title("Try out colors!")

color = st.color_picker("Pick a color", "#FF5733")
st.write(f"You selected: {color}")

st.markdown(
    f'<div style="width:100px;height:40px;background-color:{color};border-radius:8px;border:1px solid #888;"></div>',
    unsafe_allow_html=True
)

import base64
"""### gif from local file"""
file_ = open("/home/gionmattia/Downloads/animated_marta2.gif", "rb")
contents = file_.read()
data_url = base64.b64encode(contents).decode("utf-8")
file_.close()

st.markdown(
    f'<img src="data:image/gif;base64,{data_url}" width="150" alt="cat gif">',
    unsafe_allow_html=True,
)

import base64

file_ = open("/home/gionmattia/Downloads/animated_marta2.gif", "rb")
contents = file_.read()
data_url = base64.b64encode(contents).decode("utf-8")
file_.close()

st.markdown(
    f'<img src="data:image/gif;base64,{data_url}" width="150" alt="cat gif">',
    unsafe_allow_html=True,
)

st.write("""RiboMARTA - Roadmap for reference.
- IDEAS: introduce the X-squared breakdowns.
- IDEAS: introudce the effieniency method (after testing).
- Introduce a "select organism" function that loads organisms. Load one, cache, forget about it.
- ADD ORF PLOTS (so load and reconstruct sequence from FASTA)
- Implement auto-noise selection (picks a range of low z-scores from transcript).
- Option to run test in a different way, if 3'UTR not present or too small (was it wilcoxon?).
- Implement triplet periodicity in the visualisations (that would be important...)?
- At the present time, I consider the stop codon as "part of the CDS"..that's wrong. Need to fix it.
""")

st.write(
    f"utr5 {st.session_state.utr5_length}, "
    f"cds {st.session_state.cds_length}, "
    f"utr3 {st.session_state.utr3_length}"""
)