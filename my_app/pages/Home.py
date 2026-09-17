# Home.py
import streamlit.elements.widgets.select_slider

from home_helpers import *
from app_helpers import *



# Initialise runs dictionary, to store multiple runs. Each run is a AnalysisStatus class, with inputs, outputs, and logical flags.
if "analysis_runs" not in st.session_state:
    st.session_state.analysis_runs ={}


##### Welcome and Introduction messages) #####
st.title("Welcome to RiboMARTA")
col1, col2 = st.columns([1,4])
with col1:
    st.image("martas/marta.png", width = 200)
with col2:
    st.markdown("""Welcome to the first prototype for RiboMARTA. My name is Marta, and I am here to guide you
through the app and the analysis steps. You are currently in the starting page of the analysis workflow.
    
First thing first, you need to decide which organism to use among the ones currently available 

You can then proceed to compare translation across:

- Regions of a transcript (Ensembl Transcript ID mode)
- Genomic coordinates (Genomic Coordinates mode, currently in beta)

The first mode is easier for inexperienced users, while the latter allows for more fine control over the areas tested.

- You can also upload raw data previously retrieved by RiboMARTA (File Upload mode).

Once you are settled, press 'Next' to proceed to the dataset filtering page.
""" )


##### Organism Loader #####
st.write("###  Select Organism 🦠")
st.write("To run the analysis, you need to first decide which organism (and version of their references) to use.")
organism_loader_interface()
if st.session_state.get("selected_organism") is not None:
    st.write(f"**Loaded Organism:** {st.session_state.selected_organism}")


##### Analysis Mode Selector #####
df = pd.DataFrame({
  'Analysis_options': ["None Selected", "Genomic Coordinates", "Ensembl Transcript ID", "File Upload"]
})

st.write("###  Select Analysis Mode 🔬")
st.write("Choose the analysis mode you are interested in:")
option = st.selectbox(
    'Select from list',
     df['Analysis_options'])

# Renders the different options
if option == 'Genomic Coordinates':
    genomic_coordinates_selection()

elif option == 'Ensembl Transcript ID':
    ensembl_txid_selection()

elif option ==  "File Upload":
    upload_selection()

elif option ==  "None Selected":
    st.write("To move forward, please select the analysis mode you desire.")