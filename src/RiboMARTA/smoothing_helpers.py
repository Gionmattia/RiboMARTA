import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import zscore

def genco_data_smoother(current_mode, option: str = "No Smoothing"):

    if current_mode.aggregation is None:
        # Resets dataframe to original, untouched version
        current_mode.genomic_merged = current_mode.genomic_merged_backup.copy()
    elif current_mode is not None:
        # Restore to post-aggregation version, if aggregation was done
        current_mode.genomic_merged = current_mode.post_aggregation_backup.copy()

    if option == "No Smoothing":
        current_mode.smoothing = False
        return

    elif option == "Mappability Smoothing (Proportional)": # Mappability Smoothing (Proportional)

        # do prop smoothing
        df = current_mode.genomic_merged
        map_df = current_mode.genomic_merged_mappability_df

        # Save position column
        position = df["Position"]
        data_cols = df.drop(columns="Position")

        # Smooth data
        smoothed = data_cols.mul(map_df["Mappability"], axis=0)

        # Reassemble
        prop_smoothed_df =  pd.concat([position, smoothed], axis=1)

        # reassign
        current_mode.genomic_merged = prop_smoothed_df
        return

    elif option == "Mappability Smoothing (Aggressive)": #Mappability Smoothing (Aggressive)

        # do smoothing
        df = current_mode.genomic_merged
        map_df = current_mode.genomic_merged_mappability_df.copy()
        map_df.loc[map_df['Mappability'] < 1, 'Mappability'] = 0

        # Save position column
        position = df["Position"]
        data_cols = df.drop(columns="Position")

        # # Multiply only the data columns
        # smoothed = data_cols.mul(map_df["Mappability"], axis=0)
        #
        # agg_smoothed_df = pd.concat([position, smoothed], axis=1)

        # Create mask for rows where mappability < 1
        low_map = (map_df["Mappability"] < 1).to_numpy()
        data_cols.loc[low_map, :] = np.nan

        # assemble result
        agg_smoothed_df = pd.concat([position, data_cols], axis=1)
        # reassign
        current_mode.genomic_merged = agg_smoothed_df
        return

    elif option == "Smooth Criminal (Experimental)":  # Smooth Criminal (Experimental)

        z_thresh = current_mode.smoothing_criminal_thresh

        # slice the genomic merged dataframe.
        df = current_mode.genomic_merged
        if "Position" in df.columns:
            df = df.drop("Position", errors="ignore")

        x_bound = current_mode.genomic_merged_ranges["x_slice"]
        y_bound = current_mode.genomic_merged_ranges["y_slice"]
        n_bound = current_mode.genomic_merged_ranges["n_slice"]

        # Using the slices, which have already been stored individually
        df_x = df[x_bound[0]: x_bound[1]].copy()
        df_y = df[y_bound[0]: y_bound[1]].copy()
        df_n = df[n_bound[0]: n_bound[1]].copy()

        # run smooth criminal 3 times (one per slice)
        annie_x = smooth_criminal(df_x, z_thresh)
        annie_y = smooth_criminal(df_y,z_thresh)
        annie_n = smooth_criminal(df_n, z_thresh)

        # reassemble
        criminal_smoothed_df = pd.concat([annie_x, annie_y, annie_n])

        # reassign to genomic merged
        current_mode.genomic_merged = criminal_smoothed_df
        return


def genco_data_aggregator(current_run, aggregation):

    # Restores to original backup, undoing all previous modifications
    current_run.genomic_merged = current_run.genomic_merged_backup.copy()
    # Removes any flag from prior smoothing, as it's resetting the data
    current_run.smoothing = None
    current_run.smoothing_criminal_thresh = None

    # Resets everything
    if aggregation == None:
        current_run.meta_agg_option = None
        current_run.post_aggregation_backup = None # Removes post_aggregation_backup
        return

    # Metadata Aggregation
    elif aggregation  ==  "Aggregate by Metadata":

        # Consistency, remove position (though it won't be aggregated, but better safe than sorry)
        current_run.genomic_merged = current_run.genomic_merged.drop(columns = "Position")
        metadata = st.session_state.metadata.copy()
        #subset to keep only selected meta
        sel_studies = current_run.genomic_merged.columns
        metadata = metadata[metadata["Run"].isin(sel_studies)]

        mapping = metadata.set_index("Run")[current_run.meta_agg_option]

        current_run.members_per_group = mapping.value_counts().to_frame()
        current_run.members_per_group["Group"] = current_run.members_per_group.index
        current_run.members_per_group.index.name = None
        current_run.genomic_merged = current_run.genomic_merged.T.groupby(mapping).sum().T
        current_run.genomic_merged["Position"] = current_run.genomic_merged_backup["Position"]

        current_run.post_aggregation_backup = current_run.genomic_merged.copy()
        return

    # Total Aggregation
    elif aggregation  == "Total Aggregation":

        current_run.genomic_merged = current_run.genomic_merged.drop(columns = "Position")
        studies_aggregated = current_run.genomic_merged.shape[1]
        current_run.genomic_merged = current_run.genomic_merged.sum(axis=1).to_frame("Global Aggregate")

        current_run.genomic_merged["Position"] = current_run.genomic_merged_backup["Position"]

        current_run.post_aggregation_backup = current_run.genomic_merged.copy()
        current_run.members_per_group = pd.DataFrame([["Global Aggregate", studies_aggregated]], columns=["Group", "count"])
        return


# Calls the appropriate smoothing function, depending on the option given.
def txid_data_smoother(current_mode, option: str = "No Smoothing"):

    #CHECK
    if current_mode.mappability_df is not None:
        assert current_mode.all_studies_df.index.equals(current_mode.mappability_df.index), "Index mismatch between studies and mappability"

    # New bit: allows to do smoothing ON TOP of aggregation
    # RESTORE BACKUP
    if current_mode.aggregation is None:
        # Resets dataframe to original, untouched version
        current_mode.all_studies_df = current_mode.all_studies_df_backup.copy()
    elif current_mode.aggregation is not None:
        # Restore to post-aggregation version, if aggregation was done
        current_mode.all_studies_df = current_mode.post_aggregation_backup.copy()

    if option == "No Smoothing":
        current_mode.smoothing = None
        return

    elif option == "Mappability Smoothing (Proportional)": # Mappability Smoothing (Proportional)

        # do prop smoothing
        prop_smoothed_df = prop_smoother(current_mode)
        # reassign
        current_mode.all_studies_df = prop_smoothed_df
        return

    elif option == "Mappability Smoothing (Aggressive)": #Mappability Smoothing (Aggressive)

        # do something
        agg_smoothed_df = agg_smoother(current_mode)
        # reassign
        current_mode.all_studies_df = agg_smoothed_df
        return

    elif option == "Smooth Criminal (Experimental)":  # Smooth Criminal (Experimental)

        df = current_mode.all_studies_df
        z_thresh = current_mode.smoothing_criminal_thresh

        utr5_length = current_mode.utr5_length
        cds_length = current_mode.cds_length

        # If annotation info is missing (could happen for lncRNAs) then everything smoothed together
        if utr5_length is None or cds_length is None:
            st.warning(
                "UTR or CDS length not available — applying criminal smoothing to "
                "the full transcript instead of region-wise."
            )
            return smooth_criminal(df, z_thresh)


        # subset the dataframe in three regions
        df_5 = df.iloc[:utr5_length].copy()
        df_cds = df.iloc[utr5_length: utr5_length+cds_length].copy()
        df_3 = df.iloc[utr5_length+cds_length:].copy()

        # run smooth criminal 3 times (one per region)
        annie_5 = smooth_criminal(df_5, z_thresh)
        annie_cds = smooth_criminal(df_cds,z_thresh)
        annie_3 = smooth_criminal(df_3, z_thresh)

        # reassemble
        criminal_smoothed_df = pd.concat([annie_5, annie_cds, annie_3])
        assert criminal_smoothed_df.index.equals(df.index)

        # reassign
        current_mode.all_studies_df = criminal_smoothed_df
        return


def prop_smoother(current_mode):

    df = current_mode.all_studies_df
    map_df = current_mode.mappability_df
    prop_smoothed_df = df.mul(map_df["Mappability"], axis = 0)

    return prop_smoothed_df


def agg_smoother(current_mode):

    agg_smoothed_df  = current_mode.all_studies_df.copy()
    map_df = current_mode.mappability_df.copy()
    # map_df.loc[map_df['Mappability'] < 1, 'Mappability'] = 0
    # agg_smoothed_df = df.mul(map_df["Mappability"], axis=0)

    # Create mask for rows where mappability < 1
    low_map = (map_df["Mappability"] < 1).to_numpy()
    agg_smoothed_df .loc[low_map, :] = np.nan

    return agg_smoothed_df

# NDR We appply capping only to high values because, if we were to do symmetrically, sutdies with strong outliers might
# End up having 0-count positions end in the negative space of the z-score -> those values mightget replaced with higher ones, thus
# generating signal where there should be none! We do not want an erroneous iflation.
def smooth_criminal(df, z_thresh):

    # log transform and get z-score per column
    df_log2 = np.log2(df + 1)
    df_z = df_log2.apply(zscore, ddof=1)

    # Find high outliers (criminals)
    high_outlier_mask = df_z > z_thresh

    # Per-column capping. Yes, we are doing a "numerical glass ceiling"
    highest_allowed = df.where(~high_outlier_mask).max()

    # Replace
    criminal_smoothed_df = df.mask(high_outlier_mask, highest_allowed, axis=1)

    return criminal_smoothed_df


def data_aggregator(current_run, aggregation):

    # Restores to original backup, undoing all previous modifications
    current_run.all_studies_df = current_run.all_studies_df_backup.copy() # current_run.all_studies_df_pre__aggregation_backup.copy()
    # Removes any flag from prior smoothing, as it's resetting the data
    current_run.smoothing = None
    current_run.smoothing_criminal_thresh = None

    # Resets everything
    if aggregation == None:
        current_run.meta_agg_option = None
        current_run.post_aggregation_backup = None # Removes pre_aggregation_backup
        return

    # Metadata Aggregation
    elif aggregation  ==  "Aggregate by Metadata":

        # Consistency, remove position (though it won't be aggregated, but better safe than sorry)
        current_run.all_studies_df = current_run.all_studies_df.drop(columns = "Position")
        metadata = st.session_state.metadata.copy()
        #subset to keep only selected meta
        sel_studies = current_run.all_studies_df.columns
        metadata = metadata[metadata["Run"].isin(sel_studies)]

        mapping = metadata.set_index("Run")[current_run.meta_agg_option]

        current_run.members_per_group = mapping.value_counts().to_frame()
        current_run.members_per_group["Group"] = current_run.members_per_group.index
        current_run.all_studies_df = current_run.all_studies_df.T.groupby(mapping).sum().T
        current_run.all_studies_df["Position"] = current_run.all_studies_df_backup["Position"]   # current_run.all_studies_df_pre__aggregation_backup["Position"]

        current_run.post_aggregation_backup = current_run.all_studies_df.copy()
        return

    # Total Aggregation
    elif aggregation  == "Total Aggregation":

        current_run.all_studies_df = current_run.all_studies_df.drop(columns = "Position")
        studies_aggregated = current_run.all_studies_df.shape[1]
        current_run.all_studies_df = current_run.all_studies_df.sum(axis=1).to_frame("Global Aggregate")

        current_run.all_studies_df["Position"] = current_run.all_studies_df_backup["Position"]   # current_run.all_studies_df_pre__aggregation_backup["Position"]

        current_run.post_aggregation_backup = current_run.all_studies_df.copy()
        current_run.members_per_group = pd.DataFrame([["Global Aggregate", studies_aggregated]], columns=["Group", "count"])
        return