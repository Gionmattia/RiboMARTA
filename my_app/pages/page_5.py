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

ranges_selector(200, "x")

ranges_selector(200, "y")