# Home.py

from home_helpers import *
from app_helpers import *

# Initialise runs dictionary, to store multiple runs. Each run is a AnalysisStatus class, with inputs, outputs, and logical flags.
if "analysis_runs" not in st.session_state:
    st.session_state.analysis_runs ={}

st.title("Welcome to RiboMARTA")
col1, col2 = st.columns([1,4])
with col1:
    st.image("martas/marta.png", width = 200)
with col2:
    st.write("""Welcome to the first prototype for RiboMARTA. My name is Marta, and I am here to guide you
    through the app and the analysis steps.""")
    
    st.write("""If Pasha likes me, I will be allowed to live.""")
    
    st.write("""You are currently in the Home pag. First thing first, you need to decide which organism to use,
    among the ones currently available.
    You can then proceed to compare translation across regions, which can be extracted from their 
    genomic coordinates, an ensembl transcript ID, or by uploading raw data from a previous run fo RiboMARTA.
     Once you are settled, press
    'Next'.""")

## Organism Loader
st.write("###  Select Organism 🦠")
st.write("To run the analysis, you need to first decide which organism (and version of their references) to use.")
organism_loader_interface()
# st.session_state.meta_filter = st.toggle("Do you wish to filter studies by available metadata, before the analysis?")

# Sets up possible choices for analysis
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




