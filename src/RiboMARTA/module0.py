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
#import pyranges as pr
import pandas as pd
from joblib import Parallel, delayed
import os
import multiprocessing as mp
mp.set_start_method("fork", force=True)
import random
import time


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
    :return: The dataframe with transcript_ranges to extract from genomic
    """
    data = []
    region_names = ["Baseline Region", "Target Region", "Noise Region"]
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
    This is the function used to extract coverage info when the user specifies the three genomic transcript_ranges manually.

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

    pre_feature_x = final["Baseline Region"]
    post_feature_y = final["Target Region"]
    noise_n = final["Noise Region"]

    # Replace na with 0s instead
    pre_feature_x.fillna(0, inplace=True)
    post_feature_y.fillna(0, inplace=True)
    noise_n.fillna(0, inplace=True)

    return pre_feature_x, post_feature_y, noise_n


# General-use Helper
def global_aggregate(input_df: pd.DataFrame,
                     normalize: bool = False,
                     keep_pos: bool= False
                     ) -> pd.DataFrame:
    """
    Compute global aggregate coverage across all samples.

    Sums coverage values across all runs (columns) for each position (row).
    Optionally normalizes each sample by its total coverage before aggregating.

    Args:
        input_df: DataFrame with coverage data. Rows are positions, columns are
                  sample IDs (SRRs), must contain a 'Position' column
        normalize: If True, normalize each sample by total coverage before
                   aggregating; if False, sum raw counts (default: False)
        keep_pos: If True, include 'Position' column in output (default: False)

    Returns:
        DataFrame with 'Global Aggregate' column containing summed (or normalized
        and summed) coverage per position. Optionally includes 'Position' column.

    Examples:
        >>> # Raw count aggregation
        >>> ga_raw = global_aggregate(coverage_df, normalize=False)

        >>> # Normalized aggregation (each sample contributes equally)
        >>> ga_norm = global_aggregate(coverage_df, normalize=True)
    """
    # Work on a copy to avoid modifying original
    df = input_df.copy()

    # Extract position info if needed
    positions = df["Position"] if keep_pos and "Position" in df.columns else None

    # Remove Position column and handle missing values
    if "Position" in df.columns:
        df.drop(columns="Position", inplace=True)
    df.fillna(0, inplace=True)

    if normalize:
        # Normalize each sample by its total coverage
        col_sums = df.sum()
        # Replace 0 with 1 to avoid division by zero (results in 0/1 = 0)
        df = df / col_sums.replace(0, 1)


    # Computes global aggregate (ignores Sequence column)
    result_df = pd.DataFrame({
        "Global Aggregate": df.select_dtypes(include='number').sum(axis=1)
    })

    # Add position column if requested and present in the original input_df
    if keep_pos and "Position" in input_df.columns:
        result_df["Position"] = positions

    return result_df