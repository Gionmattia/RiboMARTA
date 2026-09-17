import streamlit as st

def txid_mode_style():
    st.markdown("""
    <style>
        /* 1. The Base Button Style (Targeting 1st and 2nd) */
        div[data-testid="stButton"]:nth-of-type(1) > button,
        div[data-testid="stButton"]:nth-of-type(2) > button,
        div[data-testid="stDownloadButton"] > button {
            background-color: #000000;
            color: #ffffff;
            border: 2px solid #555555;
            border-radius: 10px;
            padding: 0.5em 1em;
            transition: all 0.3s ease;
        }
    
        /* 2. The Hover State (Targeting 1st and 2nd) */
        div[data-testid="stButton"]:nth-of-type(1) > button:hover,
        div[data-testid="stButton"]:nth-of-type(2) > button:hover,
        div[data-testid="stDownloadButton"] > button:hover  {
            border-color: #FFA500;
            color: #FFA500;
            background-color: #1a1a1a;
        }
    
        /* 3. The Clicked/Active State (Targeting 1st and 2nd) */
        div[data-testid="stButton"]:nth-of-type(1) > button:active,
        div[data-testid="stButton"]:nth-of-type(2) > button:active,
        div[data-testid="stDownloadButton"] > button:active {
            background-color: #FFA500;
            color: #ffffff;
            border-color: #FFA500;
        }
    
        /* 4. Focus State (Targeting 1st and 2nd) */
        div[data-testid="stButton"]:nth-of-type(1) > button:focus:not(:active),
        div[data-testid="stButton"]:nth-of-type(2) > button:focus:not(:active),
        div[data-testid="stDownloadButton"] button:focus:not(:active) >
         {
            border-color: #FFA500;
            color: #FFA500;
            box-shadow: none;
        }
    
    
    /* ===== EXPANDERS ===== */
    /* New versions */
    div[data-testid="stExpander"] {
        background-color: #2b2b2b !important;
        border-radius: 8px;
        border: 1px solid #444;
    }
    
    div[data-testid="stExpander"] > details {
        background-color: #f0aa2d !important;
    }
    
    div[data-testid="stExpander"] > details > summary,
    div[data-testid="stExpander"] > details > summary p,
    div[data-testid="stExpander"] > details > summary span {
        background-color: #f0aa2d !important;
        color: #000000 !important;
        padding: 0.17rem 0.5rem;
        font-size: 1.15rem !important;
    }
    
    div[data-testid="stExpander"] > details[open] > summary {
        background-color: #f0aa2d !important;
    }
    
    /* Content inside expander */
    div[data-testid="stExpanderDetails"] {
        background-color: #2b2b2b !important;
        padding: 1rem;
    }
    
    /* ===== TABS ===== */
    button[data-baseweb="tab"] {
        background-color: #f0aa2d !important;
        color: #000000 !important;
        padding: 0.75rem 1.5rem !important;
        font-size: 1.15rem !important;
        font-weight: 600 !important;
    }
    button[data-baseweb="tab"] p {
        font-size: 1.1rem !important;
    }
    
    button[data-baseweb="tab"][aria-selected="true"] {
        background-color: #f0aa2d !important;
        color: #000000 !important;
    }
    
    div[data-baseweb="tab-highlight"] {
        background-color: #f0c694 !important;
        height: 0.45rem !important;
    }
    
    div[data-baseweb="tab-panel"] {
        background-color: #2b2b2b !important;
        padding: 1.5rem;
        border-radius: 0 0 8px 8px;
    }
    
    div[data-baseweb="tab-list"] {
        background-color: transparent !important;
        gap: 4px;
    }
    /* Unselected radio - default state (outer ring) */
    div[data-testid="stRadio"] label[data-baseweb="radio"] > div:first-child 
        border: 2px solid #888888 !important;
        background-color: transparent !important;
        box-sizing: border-box !important;
    }
    /* Unselected radio - hover state (outer ring) */
    div[data-testid="stRadio"] label[data-baseweb="radio"]:hover:not(:has(input:checked)) > div:first-child {
        border: 2px solid #d0a454 !important;
        background-color: rgba(30, 144, 255, 0.1) !important;
        box-sizing: border-box !important;
    }
    /* Outer ring of the selected radio */
    div[data-testid="stRadio"] label[data-baseweb="radio"]:has(input:checked) > div:first-child {
        border: 2px solid #d0a454 !important;
        background-color: transparent !important;
        box-sizing: border-box !important;
    }
    /* Inner dot of the selected radio */
    div[data-testid="stRadio"] label[data-baseweb="radio"]:has(input:checked) > div:first-child > div {
        border: 2px solid #d0a454 !important;
        background-color: #d0a454 !important;
        box-sizing: border-box !important;
    }
    </style>
    """, unsafe_allow_html=True)
    return

def genco_mode_style():
    st.markdown("""

    <style>
        [data-testid="stTable"] thead tr th {
            background-color: #5ca274 !important;
            color: #000000 !important;        /* dark text for contrast on the green */
            font-weight: 600 !important;
        }

        /* 1. The Base Button Style (Targeting 1st and 2nd) */
        div[data-testid="stButton"]:nth-of-type(1) > button,
        div[data-testid="stButton"]:nth-of-type(2) > button,


        /* 2. The Hover State (Targeting 1st and 2nd) */
        div[data-testid="stButton"]:nth-of-type(1) > button:hover  {
            border-color: #5ca274;
            color: #5ca274;
            background-color: #1a1a1a;
        }

        /* 3. The Clicked/Active State (Targeting 1st and 2nd) */
        div[data-testid="stButton"]:nth-of-type(1) > button:active,
        div[data-testid="stButton"]:nth-of-type(2) > button:active,
        div[data-testid="stDownloadButton"] > button:active {
            background-color: #78c38c;
            color: #ffffff;
            border-color: #78c38c;
        }

        /* 4. Focus State (Targeting 1st and 2nd) */
        div[data-testid="stButton"]:nth-of-type(1) > button:focus:not(:active),
        div[data-testid="stButton"]:nth-of-type(2) > button:focus:not(:active),
        div[data-testid="stDownloadButton"] button:focus:not(:active) >
         {
            border-color: #FFA500;
            color: #FFA500;
            box-shadow: none;
        }


    /* ===== EXPANDERS ===== */
    /* New versions */
    div[data-testid="stExpander"] {
        background-color: #2b2b2b !important;
        border-radius: 8px;
        border: 1px solid #444;
    }

    div[data-testid="stExpander"] > details {
        background-color: #5ca274 !important;
    }

    div[data-testid="stExpander"] > details > summary,
    div[data-testid="stExpander"] > details > summary p,
    div[data-testid="stExpander"] > details > summary span {
        background-color: #5ca274 !important;
        color: #000000 !important;
        padding: 0.17rem 0.5rem;
        font-size: 1.15rem !important;
    }

    div[data-testid="stExpander"] > details[open] > summary {
        background-color: #5ca274 !important;
    }

    /* Content inside expander */
    div[data-testid="stExpanderDetails"] {
        background-color: #2b2b2b !important;
        padding: 1rem;
    }

    /* ===== TABS ===== */
    button[data-baseweb="tab"] {
        background-color: #5ca274 !important;
        color: #000000 !important;
        padding: 0.75rem 1.5rem !important;
        font-size: 1.15rem !important;
        font-weight: 600 !important;
    }
    button[data-baseweb="tab"] p {
        font-size: 1.1rem !important;
    }

    button[data-baseweb="tab"][aria-selected="true"] {
        background-color: #5ca274 !important;
        color: #000000 !important;
    }

    div[data-baseweb="tab-highlight"] {
        background-color: #78c38c !important;
        height: 0.45rem !important;
    }

    div[data-baseweb="tab-panel"] {
        background-color: #2b2b2b !important;
        padding: 1.5rem;
        border-radius: 0 0 8px 8px;
    }

    div[data-baseweb="tab-list"] {
        background-color: transparent !important;
        gap: 4px;
    }
    /* Unselected radio - default state (outer ring) */
    div[data-testid="stRadio"] label[data-baseweb="radio"] > div:first-child {
        border: 2px solid #888888 !important;
        background-color: transparent !important;
        box-sizing: border-box !important;
    }
    /* Unselected radio - hover state (outer ring) */
    div[data-testid="stRadio"] label[data-baseweb="radio"]:hover:not(:has(input:checked)) > div:first-child {
        border: 2px solid #00855d !important;
        background-color: rgba(30, 144, 255, 0.1) !important;
        box-sizing: border-box !important;
    }
    /* Outer ring of the selected radio */
    div[data-testid="stRadio"] label[data-baseweb="radio"]:has(input:checked) > div:first-child {
        border: 2px solid #00855d !important;
        background-color: transparent !important;
        box-sizing: border-box !important;
    }
    /* Inner dot of the selected radio */
    div[data-testid="stRadio"] label[data-baseweb="radio"]:has(input:checked) > div:first-child > div {
        border: 2px solid #00855d !important;
        background-color: #00855d !important;
        box-sizing: border-box !important;
    }
    </style>
    """, unsafe_allow_html=True)

# Style bubble - Used to recap inputs/options throughout the app
def bubble(text, bg_color="#000000", outline_color="#f0aa2d", text_color="#f0aa2d", outline_width=1):
    st.markdown(
        f'<div style="'
        f'background:{bg_color};'
        f'color:{text_color};'
        f'border:{outline_width}px solid {outline_color};'
        f'padding:10px 14px;'
        f'border-radius:10px;'
        f'display:inline-block;'
        f'font-weight:600;'
        f'">{text}</div>',
        unsafe_allow_html=True,
    )
