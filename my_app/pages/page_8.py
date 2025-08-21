import streamlit as st
from transcriptID_mode_helpers import *
import yaml
from app_helpers import *

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

current_run = st.session_state.runs[st.session_state.quirky_id]
st.write(f"######  :gray[ Analysis ID: {st.session_state.quirky_id}]")
st.write(st.session_state.runs.keys())


# st.write("""RiboMARTA - Roadmap for reference.
# - IDEAS: introduce the X-squared breakdowns.
# - IDEAS: introudce the effieniency method (after testing).
# - Implement auto-noise selection (picks a range of low z-scores from transcript).
# - Option to run test in a different way, if 3'UTR not present or too small (was it wilcoxon?).
# """)

if st.button("Plese, report!"):
    info = compile_report(current_run)
    st.write(info)



    # Your multi-line text, e.g., a report of parameters
    report_text = """
    Report Date: August 21, 2025
    -----------------------------
    Model: XYZ-123
    Accuracy: 95.2%
    Precision: 94.8%
    Recall: 96.1%
    F1-Score: 95.4%
    Parameters:
      - Learning Rate: 0.01
      - Batch Size: 32
      - Epochs: 50
    """

    st.text_area(
        "Parameter Report",
        report_text,
        height=250
    )