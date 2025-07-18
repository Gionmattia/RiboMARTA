import streamlit as st
from transcriptID_mode_helpers import *
import yaml

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

st.write("## Please, select the organism to work with")
with open("config.yaml") as f:
    configs = yaml.safe_load(f)

st.write(list(configs.keys()))


# st.write("""RiboMARTA - Roadmap for reference.
# - IDEAS: introduce the X-squared breakdowns.
# - IDEAS: introudce the effieniency method (after testing).
# - Introduce a "select organism" function that loads organisms. Load one, cache, forget about it.
# - Implement auto-noise selection (picks a range of low z-scores from transcript).
# - Option to run test in a different way, if 3'UTR not present or too small (was it wilcoxon?).
# """)
