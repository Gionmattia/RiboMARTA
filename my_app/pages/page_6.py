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

coordinates =   {"x_slice": (1,10),
                "y_slice": (1,10),
                "n_slice": (1,10)}
st.write(f"""Reminder! You selected the following ranges for the analysis:

- First Region - From {coordinates["x_slice"][0]} to {coordinates["x_slice"][1]}
- Second Region - From {coordinates["y_slice"][0]} to {coordinates["y_slice"][1]}
- Noise control Region - From {coordinates["n_slice"][0]} to {coordinates["n_slice"][1]}
""")

st.session_state.my_ranges[0]