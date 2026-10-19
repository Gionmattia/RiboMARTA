## RiboMARTA - <u>Ribo</u>-Seq  <u>m</u>RNA <u>A</u>nalysis of <u>R</u>egional <u>T</u>ranslational <u>A</u>ctivity

Welcome to RiboMARTA. This app provides an intuitive interface for the comparative analysis of regional translational 
activity, using Ribosome profiling (Ribo-Seq) data.


## Getting Started

### Conda

You can create an empty conda environment and then install RiboMARTA from GitHub using the command


```bash
conda create -n ribomarta python=3.10
conda activate ribomarta
pip install git+https://github.com/Gionmattia/RiboMARTA.git
```

### Input files required

RiboMARTA requires the following files in order to run, which are all passed through a single .yml/.yaml config file:

NB. The file must contain the **absolute paths** to the files

| Input | Format | Required | Notes |
|---|---|---|---|
| Genome sequence | FASTA (`.fa`) | Yes | Same version used for the generation of the BigWigs.|
| Genome annotation | Ensembl GTF (`.gtf`) | Yes | Same version used for the generation of the BigWigs.|
| Ribo-Seq data | Folder of BigWig (`.bw`) files | Yes | One file per dataset, containing inferred A- or P-site positions. Each file must be named after its SRA Run accession (e.g. `SRR1234567.bw`). |
| Metadata | CSV (`.csv`) | Yes | Must include a `Run` column with the SRR accession of every dataset. All other fields are optional. |
| Mappability track | UCSC Umap k24 as a BigWig (`.bw`) | No | Enables display of a mappability track and its use during analysis. |

### Notes on metadata

The only enforced requirement is that every dataset is identified by its SRR accession in the `Run` column. 
More accurate metadata (e.g. cell type, treatment, condition) is strongly recommended, as it makes results easier to
interpret. Runs with a missing value in a given field are excluded from any analysis that uses that field.

### Running

You can run RiboMARTA using the command 

```bash
conda activate ribomarta
ribomarta --config /absolute/path/to/config.yaml
```

If no config file is given, RiboMARTA will default to the built-in example config.

You can run RiboMARTA through streamlit too (not sure why you'd like to, but anyway):

```bash
conda activate ribomarta
streamlit run <path/to/RiboMARTTA_main.py>  -- --config /absolute/path/to/config.yaml
```
Again, if no .yml is given, it will default to the built-in example.


## Cite 

Author: Gionmattia Carancini

DOI: 10.5281/zenodo.23477139

