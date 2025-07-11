"""
module0.py

This module contains all the functions used to extract the coverage values from a (genomic) bigwig.
It contains the functions to extract three uninterruped regions from the coordinates, as well as
the functions to extract a single transcript.
and to calculate the ratios.
Author: Gionmattia Carancini
Date: 26-05-2025
"""

# Import libraries
import pyBigWig
import pyranges as pr
import pandas as pd
from joblib import Parallel, delayed
import os
import multiprocessing as mp
mp.set_start_method("fork", force=True)
import plotly.graph_objects as go
import random
import re
import time

# # OLD VERSION - BACKUP Helper - Extracts the coverage info from a single bigwig
# def _extract_one(bw_path, regions_df):
#     """
#     :param bw_path: path to a signle .bigwig file
#     :param regions_df: dataframe with genomic regions to extract
#     :return: tuple. Study ID and dictionary of the three regions dataframes.
#     """
#
#     bw   = pyBigWig.open(bw_path)
#     sid  = os.path.basename(bw_path).split(".")[0]
#     out  = {}
#
#     # Convert the coordinates from 1-indexed inclusive to BigWig's 0-indexed end-exclusive [,)
#     # Not necessary, pyranges already converts the coordinates system
#     # regions_df["Start"] = regions_df["Start"] - 1
#
#     # iterate directly over the DataFrame
#     for row in regions_df.itertuples(index=False):
#         chrom, start, end, region = (
#             row.Chromosome, int(row.Start), int(row.End), row.Region
#         )
#         vals = bw.values(chrom, start, end, numpy=True)
#         out[region] = pd.DataFrame({
#             "Position": range(start, end),
#             sid: vals
#         })
#
#     bw.close()
#     return sid, out

# Helper - Extracts the coverage info from a single bigwig. Used in both Genomic Coordinates and TranscriptID modes
def _extract_one(bw_path, regions_df):
    """
    :param bw_path: path to a signle .bigwig file
    :param regions_df: dataframe with genomic regions to extract
    :return: tuple. Study ID and dictionary of the three regions dataframes.
    """

    sid  = os.path.basename(bw_path).split(".")[0]
    out  = {}

    attempts = 0
    max_attempts = 2

    while attempts < max_attempts:
        try:
            bw = pyBigWig.open(bw_path)
            # iterate directly over the DataFrame
            for row in regions_df.itertuples(index=False):
                chrom, start, end, region = (
                    row.Chromosome, int(row.Start), int(row.End), row.Region
                )
                vals = bw.values(chrom, start, end, numpy=True)
                out[region] = pd.DataFrame({
                    "Position": range(start, end),
                    sid: vals
                })
            bw.close()
            break  # success, so exit the loop

        # Failsafe. Sometimes parallelisation cannot work if same file is opened at the same time by != forks.
        # We introduce a random sleep time for a worker, just to resolve the issue. If that does not work, then....
        # returns 0s instead.
        except:
            # Increase attempts counter
            attempts += 1
            if attempts < max_attempts:
                # Wait a random time to ensure no overlap
                wait = random.uniform(0.3, 1.5)
                time.sleep(wait)
            # If issue cannot be resolved by waiting, returns zeros.
            else:
                for row in regions_df.itertuples(index=False):
                    chrom, start, end, region = (
                        row.Chromosome, int(row.Start), int(row.End), row.Region
                    )
                    # Fill with zeros (or None if you prefer)
                    out[region] = pd.DataFrame({
                        "Position": range(start, end),
                        sid: None
                    })

    return sid, out


# Function to format the regions_df in a standardised way (for Genomic Coordinates Mode)
def format_input_df(my_list):
    """
    Unpacks genomic coordinates and formats them into a table that is used by the code to extract coverage

    :param my_list: A list of strings (genomic coordinates) as "chrx:56,741,761-56,761,390"
    :return: The dataframe with ranges to extract from genomic
    """
    data = []
    region_names = ["First Region", "Second Region", "Noise Region"]
    for i, reg in zip(my_list, region_names):
        fri = i.replace(",", "").replace(".", "").replace(":", "-").split("-")
        chrom = fri[0]
        if chrom.lower().startswith('chr'):
            chrom = chrom[3:]
        start = int(fri[1]) - 1  # Converts from 1-indexed closed to 0-index open
        end = int(fri[2])
        data.append({"Region": reg, "Chromosome": chrom, "Start": start, "End": end})
    return pd.DataFrame(data)


# Function to extract the coverage from all bigwigs and compile the three dataframes (used in Genomic Coordinates mode).
def extract_regions_from_bw_parallel(bigwig_files, regions_df, n_jobs=1):
    """
    This is the function used to extract coverage info when the user specifies the three genomic ranges manually.

    :param bigwig_files: list of paths to .bigwig files
    :param regions_df: dataframe with genomic regions to extract
    :param n_jobs: how many cores to allocate
    :return: three dataframes, with coverage from x, y, and n regions
    """

    # Extracts in parallel the coverages from each bigwig
    # SAME AS BELOW< NEEDS TO BE FIXED. BIGWIGS SHOULD BE GIVEN IN CHUNKS TO EACH WORKER
    results = Parallel(n_jobs=n_jobs,
                       backend="multiprocessing"
                       )(
        delayed(_extract_one)(bw, regions_df) for bw in bigwig_files
    )

    # Unpacks the results and compiles everything in a single dictionary
    final = {}
    for sample_id, bw_dict in results:
        for region, df_tmp in bw_dict.items():
            if region in final:
                final[region] = final[region].merge(
                    df_tmp, on="Position", how="outer"
                )
            else:
                final[region] = df_tmp

    pre_feature_x = final["First Region"]
    post_feature_y = final["Second Region"]
    noise_n = final["Noise Region"]

    return pre_feature_x, post_feature_y, noise_n


# General-use Helper
def global_aggregate(input_df: pd.DataFrame,
                     keep_pos: bool= False
                     ) -> pd.DataFrame:
    '''
    Args:
        input_df: dataframe with coverage info per position. Rows are positions, columns are SRR Runs.
        keep_pos: whether to keep position numbers (genomic coordinate numbers) or not.

    Returns:
        pd.DataFrame with aggregated coverage per position

    Creates a new dataframe with the "Global Aggregate" column, which is just the sum of counts per position
    '''

    df =  input_df.copy()
    df.drop(columns="Position", inplace=True)
    df.fillna(0, inplace=True)

    col_sums = df.sum()
    # Replace 0 with 1s, so to not get NaN when diving by 0
    ga_df = df / col_sums.replace(0, 1)
    ga_df["Global Aggregate"] = ga_df.sum(axis = 1)
    if keep_pos == True:
        ga_df["Position"] = input_df["Position"]

    return ga_df

def genomic_mode_formatter(pre_feature_x:pd.DataFrame,
                           post_feature_y:pd.DataFrame,
                           noise_n:pd.DataFrame
                           ) -> tuple[pd.DataFrame, dict[str, tuple[int, int]]]:
    '''
    Args:
        pre_feature_x: pd.DataFrame with coverages from genomic coordinates (rows are positions, columns are SRR Runs)
        post_feature_y: pd.DataFrame with coverages from genomic coordinates (rows are positions, columns are SRR Runs)
        noise_n: pd.DataFrame with coverages from genomic coordinates
        used as noise control (rows are positions, columns are SRR Runs)

    Returns:
        genomic_merged: pd.DataFrame, the three inputs, concatenated on top of each other.
        genomic_slices: dictionary, the coordinates to re-extract the original three DataFrames.
    '''

    # Create a single df for run_analysis()
    genomic_merged = pd.concat([pre_feature_x, post_feature_y, noise_n],
                               axis=0,
                               ignore_index=True)

    length_x = len(pre_feature_x)
    length_y = len(post_feature_y)
    length_n = len(noise_n)

    # Set up the dictionary of coordinates for run_analysis()
    genomic_slices = {"x_slice": (0, length_x),
                      "y_slice": (length_x,length_x + length_y),
                      "n_slice": (length_x + length_y, length_x + length_y + length_n)}

    return genomic_merged, genomic_slices