# Home.py

from home_helpers import *
from app_helpers import *

# Initialise runs dictionary, to store multiple runs. Each run is a AnalysisStatus class, with inputs, outputs, and logical flags.
if "runs" not in st.session_state:
    st.session_state.runs ={}

st.title("Welcome to RiboMARTA")
col1, col2 = st.columns([1,4])
with col1:
    st.image("martas/marta.png", width = 200)
with col2:
    st.write("""Welcome to the first draft of RiboMARTA. My name is Marta, and I am here to guide you
    through the steps of the analysis.""")
    
    st.write("""If Pasha likes me, I will be allowed to live.""")
    
    st.write("""You are currently in the Home page of the analysis. In here you can decide whether to extract regions using
    genomic coordinates, an ensembl transcript ID, or use raw data from a previous run fo RiboMARTA. Once you are settled, press
    'Next'.""")

## Organism Loader (wrapped up as a function)
st.write("###  Select Organism 🦠")
st.write("To run the analysis, you need to first decide which organism (and version of their references) to use.")
organism_loader_interface()

# Sets up possible choices for analysis
df = pd.DataFrame({
  'Analysis_options': ["None Selected", "Genomic Coordinates", "Ensembl Transcript ID", "File Upload"]
})

st.write("Now, choose the analysis mode you are interested in:")
option = st.selectbox(
    'Select the Analysis Mode',
     df['Analysis_options'])

# Renders the different options
if option == 'Genomic Coordinates':
    genomic_coordinates_selection()

elif option == 'Ensembl Transcript ID':
    ensembl_txid_selection()

elif option ==  "File Upload":
    upload_selection()

elif option ==  "None Selected":
    st.write("Please select the analysis mode you desire.")
    # session_state clean-up (in this case, removes everything)
    # keys_to_reset = ["my_ranges", "picture", "txid", "genomic_slices", "all_studies_df", "utr5_length",
    #                          "cds_length", "utr3_length"]
    # reset_session_state(keys_to_reset)



