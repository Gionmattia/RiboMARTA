# transcriptID_mode_helpers.py

import streamlit as st

def ranges_selector(full_mRNA_length, region):
    """
    A simple function to automatically render the selectors for ranges
    :param full_mRNA_length: Intuitively, the full cDNA length
    :param region: internal, needs to be different every time this function is called, or input keys will overlap and cause error
    :return: a tuple, with two int defining a range on the mRNA
    """
    x_start = st.number_input("Transcript region start coordinate (es: 10)",
                              min_value=0,
                              max_value=full_mRNA_length,
                              step=1,
                              key = f"start_{region}",
                              )
    x_end = st.number_input("Transcript region end coordinate (es: 100)",
                            min_value=0,
                            max_value=full_mRNA_length,
                            step=1,
                            key=f"end_{region}",
                            )
    st.write(f"Selected region {x_start} to {x_end}")

    return (x_start, x_end)

    # st.session_state.x_slice = (x_start, x_end)

