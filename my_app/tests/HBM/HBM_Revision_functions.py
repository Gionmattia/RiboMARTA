# This bit used to create the input data. I will wrap it as a function.
import pandas as pd
import pymc as pm
import numpy as np
import matplotlib.pyplot as plt
import arviz as az
import os
# Need to import preparing functions from MARTA!

#### FUNCTIONS TO REMOVE 'CAUSE IMPORTED WITH MARTA PACKAGE
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed

# needed to remove samples where we do not have info
def _clean_data(df: pd.DataFrame, verbose: bool = False) -> pd.DataFrame:
    """
    Remove runs (samples) with zero coverage across all positions.

    Filters out columns where every position has zero counts, indicating
    no sequencing data or failed library preparation.

    Args:
        df: Coverage dataframe with positions as rows, runs (SRRs) as columns
        verbose: If True, print filtering statistics (default: False)

    Returns:
        Cleaned dataframe with all-zero columns removed
    """

    n_runs_original = df.shape[1]
    # Keeps columns with at least a non-0 value.
    df_clean = df.loc[:, (df != 0).any(axis=0)]

    n_runs_after = df_clean.shape[1]

    # Originally, we used to remove all columns were x was 0...now we are doing the permutation test

    if verbose:
        n_removed = n_runs_original - n_runs_after
        print(f"Runs originally: {n_runs_original}")
        print(f"Runs with zero coverage removed: {n_removed}")
        print(f"Runs remaining: {n_runs_after}")

    return df_clean


def _permute_values(
        base_array: np.ndarray,
        noise_array: np.ndarray,
        N: int,
        random_state: int | None = None
) -> np.ndarray:
    """
    Generate null distribution by permuting baseline and noise regions.

    Shuffles the combined baseline + noise arrays N times, recomputing
    the absolute difference between means for each permutation.

    Args:
        base_array: Array of shape (n_positions_base, n_runs) from baseline region
        noise_array: Array of shape (n_positions_noise, n_runs) from noise region
        N: Number of permutations to perform
        random_state: Seed for reproducibility (default: None)

    Returns:
        Array of shape (N, n_runs) containing permuted absolute differences
    """

    # Combine the two arrays
    combined = np.concatenate((base_array, noise_array), axis=0)
    # Retrieve the original size of x
    base_size = base_array.shape[0]
    # Set the random state for permutations
    rng = np.random.default_rng(random_state)

    # Set the empty array for the results
    perms_array = np.empty((N, base_array.shape[1]))

    for i in range(N):
        # Shuffles the labels
        shuffled = np.apply_along_axis(rng.permutation, 0, combined)
        # Gets permutated sets averages
        perm_base = shuffled[:base_size].mean(axis=0)
        perm_noise = shuffled[base_size:].mean(axis=0)
        # Computes the diff between perm x and n. Stores it in the output array.
        perms_array[i, :] = abs(perm_base - perm_noise)

    return perms_array


# Computes the p-value and confidence boundary
def test_x_region(
        diffs_df: pd.DataFrame,
        perms_array: np.ndarray,
        ci: int = 95
) -> tuple[pd.DataFrame, list[str]]:
    """
    Compute p-values and identify significant runs via permutation test.

    Calculates one-tailed p-values by comparing observed differences to the
    permutation-based null distribution, then applies FDR correction.

    Args:
        diffs_df: DataFrame with observed differences, must contain 'abs_xn_diff' column
        perms_array: Array of shape (N, n_runs) with null distribution
        ci: Confidence level for FDR correction (default: 95)

    Returns:
        results: DataFrame with columns ['x avg', 'n avg', 'padj_bh', 'reject']
        good_runs: List of run IDs where reject == True
    """

    # """
    # One-tailed p-values and upper confidence bound for each series:
    #   p = (count(perms_array >= observed) + 1) / (N + 1)
    #   ci_upper = percentile(perms_array, ci)
    # """

    N = perms_array.shape[0]
    observed = diffs_df["abs_xn_diff"]

    # Now we observe how many times the perms_array abs(x-n) were more extreme than the observed values.
    greater = (perms_array >= observed.values).sum(axis=0)
    pvals = (greater + 1) / (N + 1)  # Changed from  (greater + 0.5) / (N + 1)

    # FDR correction for multiple testing
    reject, pvals_bh, _, _ = multipletests(pvals, alpha=1 - (ci / 100), method='fdr_bh')

    # Compile results
    results = pd.DataFrame({
        "x avg": diffs_df["x avg"],
        "n avg": diffs_df["n avg"],
        'padj_bh': pvals_bh,
        'reject': reject,
    }, index=observed.index)

    good_runs = results[results["reject"] == True].index.to_list()

    return results, good_runs


def test_baseline_vs_noise(
        dataframe: pd.DataFrame,
        coordinates: dict[str, tuple[int, int]],
        N: int = 10000,
        random_state: int | None = None,
        ci: int = 95
) -> tuple[pd.DataFrame, list[str]]:
    """
    Test if baseline region (x) differs significantly from noise region (n).

    Uses permutation testing to identify runs (samples) where the baseline
    region shows average raw counts significantly different from the noise region.

    Biological reasoning: not everything is espressed.
    This test makes sure that baseline is indeed expressed to a significant level.

    Args:
        dataframe: Coverage data with positions as rows, runs (SRRs) as columns
        coordinates: Dictionary with 'x_slice' and 'n_slice' keys, each containing
                    (start, end) tuples defining region boundaries
        N: Number of permutations for null distribution (default: 1000)
        random_state: Seed for reproducibility (default: None)
        ci: Confidence level for statistical testing (default: 95)

    Returns:
        results: DataFrame with columns ['x avg', 'n avg', 'padj_bh', 'reject']
                indexed by run IDs
        good_runs: List of run IDs that passed the significance test
    """

    # computes an array of observed differences between baseline(x) and noise(n).
    diffs_df, reg_x, reg_n = _compute_xn_diffs(dataframe, coordinates)

    # permutes the values between x and n, to generate Null hypothesis distribution.
    perms_array = _permute_values(reg_x, reg_n, N, random_state=random_state)

    # does the actual test, computing the p-values and adjusting them
    results, good_runs = test_x_region(diffs_df, perms_array, ci=ci)

    return results, good_runs


def _compute_xn_diffs(df: pd.DataFrame,
                      coordinates: dict[str, tuple[int, int]]
                      ) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """
    Extract baseline (x) and noise (n) regions from the dataframe, and compute their differences in terms of average raw counts.

    Args:
        df: Coverage dataframe with positions as rows, runs as columns
        coordinates: Dictionary containing 'x_slice' and 'n_slice' keys with
                    (start, end) tuples

    Returns:
        diffs_df: DataFrame with columns ['x avg', 'n avg', 'abs_xn_diff'],
                 indexed by run IDs
        reg_x: Array of shape (n_positions_x, n_runs) with x region coverage
        reg_n: Array of shape (n_positions_n, n_runs) with n region coverage
    """

    # slice out each region
    reg_x = df.iloc[coordinates["x_slice"][0]: coordinates["x_slice"][1]]
    reg_n = df.iloc[coordinates["n_slice"][0]: coordinates["n_slice"][1]]

    # Calculated as it follows so that the Pandas series keeps the rownames
    x_avg = reg_x.mean(axis=0)
    n_avg = reg_n.mean(axis=0)

    # Compute the absolute difference between x and n
    abs_xn_dff = abs(x_avg - n_avg)

    # Stores the values into a dataframe

    diffs_df = pd.DataFrame({
        "x avg": x_avg,
        "n avg": n_avg,
        "abs_xn_diff": abs_xn_dff}
    )

    return diffs_df, reg_x.to_numpy(), reg_n.to_numpy()

#### END FUNCTIONS TO REMOVE 'CAUSE IMPORTED BY MARTA ####

def prepare_input(df, coordinates):
    # Creates slices of the regions of interest: baseline (x), target (y), noise (n).
    x_slice = df.iloc[coordinates["x_slice"][0]:coordinates["x_slice"][1]].copy()
    y_slice = df.iloc[coordinates["y_slice"][0]:coordinates["y_slice"][1]].copy()
    n_slice = df.iloc[coordinates["n_slice"][0]:coordinates["n_slice"][1]].copy()

    # Generates series with Study-wise (column) sum of footprints from each slice.
    x_tot = x_slice.sum(axis=0)
    y_tot = y_slice.sum(axis=0)
    n_tot = n_slice.sum(axis=0)

    # Retrieves length of each slice:
    L_x = x_slice.shape[0]  # length baseline
    L_y = y_slice.shape[0]  # length target
    L_n = n_slice.shape[0]  # length noise

    # Sets up the dataframe with SRR - Totals per slice and their normalised version
    data = pd.DataFrame({
        "Run": x_slice.columns.values,
        'initiating': x_tot.values,  # Tot ribosomes from baseline (x)
        'finishing': y_tot.values,  # Tot ribosomes from target (y)
        'noise': n_tot.values
    })

    # For comparisons, multiply tots slice by length
    # (X/L_x and Y/L_y is equivalent to X*L_y and Y*L_x in this scenario)
    data["norm_init"] = data["initiating"] * L_y
    data["norm_finish"] = data["finishing"] * L_x

    # Stores lengths in a dict
    lengths = {"baseline": L_x, "target": L_y, "noise": L_n}

    return data, lengths


# Adds Metadata
def merge_meta(data, meta, factor):
    # Check
    if not isinstance(factor, str):
        raise TypeError(f" 'factor' must be a string, but got {type(factor).__name__} instead.")

    # merge data and metadata
    fields = ["Run", factor]
    merged_df = data.merge(
        meta[fields],
        left_on="Run",
        right_on="Run",
        how="left"
    )

    # Fill missing values (filtered later)
    for fie in fields:
        merged_df[fie] = merged_df[fie].fillna("None")

    return merged_df


# Filter Runs where the counts in target are greater than baseline for manual inspection.
def first_filter(df):
    """
    Filters based on tot counts per region. Removes Runs where norm_finish > norm_init.
    Stores those instances on df_scraped, for manual inspection.
    """

    # Creates the mask
    valid_mask = (df["norm_init"] >= df["norm_finish"])
    df_filtered = df[valid_mask]

    # Retrieves filtered_out cases
    non_valid_mask = (df["norm_init"] < df["norm_finish"])
    df_scraped = df[non_valid_mask]

    return (df_filtered, df_scraped)


# Filter Runs belonging to an analysis group/level with too little elements. We want 10 Runs per group by default.
def second_filter(df, factor, min_count=10):
    """
    Filters based on number of datasets per analysis group (factor).
    Each group/level should have at least 10 members.

    Return list of instances where a level was beneath the minimum
    """

    # Creates the mask
    too_low_mask = df[factor].value_counts() < min_count

    # Retrieves levels with too little members
    too_low = too_low_mask[too_low_mask == True].index.tolist()
    # Adds to the list missing values (Encoded previously as "None" in merge_meta()
    too_low.append("None")

    # Removes them from the data
    df_filtered = df[~df[factor].isin(too_low)].reset_index(drop=True)

    return (df_filtered, too_low)


# Compile a dictionary with the data in a format the model understand (arrays) and other mappings needed.
def format_data(df, factor):
    # prepare arrays for model
    X = df["initiating"].to_numpy()
    Y = df["finishing"].to_numpy()
    N = df["noise"].to_numpy()

    # prepare Run indexing for model (PyMC uses numbers, not strings)
    run_names = df["Run"].unique()
    run_to_idx = {name: idx for idx, name in enumerate(run_names)}
    run_idx = df["Run"].map(run_to_idx).values

    # Prepare levels indexing for model
    levels = df[factor].unique()
    levels_to_idx = {name: idx for idx, name in enumerate(levels)}
    levels_idx = df[factor].map(levels_to_idx).values

    # Prepare array where run id is matched to level id.
    run_levels_idx = (df.groupby("Run")[factor].first()
                      .reindex(run_names)
                      .map(levels_to_idx)
                      .values)  # length = n_runs

    # Store everything into a dictionary for simplicity -Will help interpreting results later
    inputs_dictionary = {
        "X": X,
        "Y": Y,
        "N": N,
        "factor": factor,
        "run_names": run_names,
        "run_to_idx": run_to_idx,
        "run_idx": run_idx,
        "levels": levels,
        "levels_to_idx": levels_to_idx,
        "levels_idx": levels_idx,
        "run_levels_idx": run_levels_idx,
        "n_obs": len(Y),
        "n_runs": len(run_names),
        "n_levels": len(levels)
    }

    return inputs_dictionary


def assemble_model(inputs_dict, lengths, belief="standard"):
    # Checks "belief" option:
    valid_beliefs = {"standard", "rare", "elevated", "agnostic"}
    if belief not in valid_beliefs:
        raise ValueError(f"belief must be one of {valid_beliefs}, got '{belief}'")

    n_levels = inputs_dict["n_levels"]
    n_runs = inputs_dict["n_runs"]
    run_idx = inputs_dict["run_idx"]
    run_levels_idx = inputs_dict["run_levels_idx"]
    levels_idx = inputs_dict["levels_idx"]
    X = inputs_dict["X"]
    L_x = lengths["baseline"]
    Y = inputs_dict["Y"]
    L_y = lengths["target"]
    N = inputs_dict["N"]
    L_n = lengths["noise"]

    with pm.Model() as model:

        ## LAYER 2: "Level" parameters (for chosen FACTOR, levels have parameters).

        # Model efficiency on logit scale
        if belief == "standard":  # Efficiencies between 1-30% (median 10%, CI 1.5-44%), like frameshifting events
            mu_factor_logit = pm.Normal('mu_factor_logit',
                                        mu=-2.2,
                                        sigma=1.0,
                                        shape=n_levels)

        elif belief == "rare":  # Efficiencies between 0.01-5% (median 1%, CI 0.05-17%), like readthrough events
            mu_factor_logit = pm.Normal('mu_factor_logit',
                                        mu=-4.6,
                                        sigma=1.5,
                                        shape=n_levels)
        elif belief == "elevated":  # Efficiencies typically 1-40% (median 10%, CI 1-53%), like selcys insertions,
            mu_factor_logit = pm.Normal('mu_factor_logit',
                                        mu=-2.2,
                                        sigma=1.2,
                                        shape=n_levels)

        elif belief == "agnostic":  # Efficiencies typically 1-40% (median 10%, CI 1-53%), like selcys insertions,
            mu_factor_logit = pm.Normal('mu_factor_logit',
                                        mu=-3.5,
                                        sigma=2.0,
                                        shape=n_levels)

        sigma_factor = pm.HalfNormal('sigma_factor',
                                     sigma=1,  # run-level variability within level
                                     shape=n_levels)

        ## LEVEL 1: Run-level efficiency (logit-normal)
        # Non-centered parameterization - models assumptions are the same, just using algebra to evade sampler issues (math magic).
        theta_run_logit_offset = pm.Normal('theta_run_logit_offset',
                                           mu=0,
                                           sigma=1,
                                           shape=n_runs)

        theta_run_logit = pm.Deterministic('theta_run_logit',
                                           mu_factor_logit[run_levels_idx] + theta_run_logit_offset * sigma_factor[
                                               run_levels_idx]
                                           )

        # Transform to (0,1) scale. Makes it interpretable again, it's a simple transformation
        theta_run = pm.Deterministic('theta_run',
                                     pm.math.invlogit(theta_run_logit))

        ## LEVEL 0: Likelihood

        # We account for noise.
        background_rate = N / L_n

        # Subtract noise from X first
        # DOUBLE CHECK pm.math.maximum with 1 flooring. Maybe is too much? Maybe 0 works best?
        clean_X = pm.math.maximum(X - background_rate * L_x, 1)  # Unnecessary, I am ready pre-filtering these cases.
        signal = theta_run[run_idx] * clean_X * (L_y / L_x)
        noise_in_Y = background_rate * L_y
        mu = signal + noise_in_Y

        phi_level = pm.Gamma("phi_level", alpha=10, beta=0.5, shape=n_levels)
        phi = phi_level[levels_idx]  # A different phi per each level of FACTOR?

        Y_obs = pm.NegativeBinomial("Y_obs", mu=mu, alpha=phi, observed=Y)

    return model


def sample_posterior(model, draws=2000, tune=1000, chains=4, target_accept=0.95, return_inferencedata=True,
                     random_seed=42):
    # Sample from posterior
    with model:
        trace_negbin = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            target_accept=target_accept,
            return_inferencedata=return_inferencedata,
            random_seed=random_seed
        )
        return trace_negbin


# Checks the summary | flags potential issues with R-hat, ESS-bulk/ESS-tail.
def check_mcmc_diagnostics(trace_negbin, n_chains=4, rhat_threshold=1.01, ess_min=None, outpath = None):
    """
    Check MCMC convergence diagnostics from an ArviZ InferenceData object.

    Parameters
    ----------
    trace_negbin : az.InferenceData
    n_chains : int
        Number of chains used (default 4). Used to compute ESS threshold.
    rhat_threshold : float
        Maximum acceptable R-hat (default 1.01, per Vehtari et al. 2021).
    ess_min : int or None
        Minimum acceptable ESS. If None, defaults to 100 * n_chains (as suggested by Stan).

    Returns
    -------
    dict with keys 'summary', 'rhat_violations', 'ess_bulk_violations',
                    'ess_tail_violations', 'passed'
    """
    if ess_min is None:
        ess_min = 100 * n_chains

    summary = az.summary(trace_negbin)

    rhat_bad = summary[summary["r_hat"] >= rhat_threshold]
    ess_bulk_bad = summary[summary["ess_bulk"] < ess_min]
    ess_tail_bad = summary[summary["ess_tail"] < ess_min]

    passed = all(df.empty for df in [rhat_bad, ess_bulk_bad, ess_tail_bad])

    # Assemble report lines
    lines = []
    lines.append(f"=== MCMC Diagnostics (thresholds: R̂ < {rhat_threshold}, ESS ≥ {ess_min}) ===")
    if passed:
        lines.append("✓ All diagnostics passed.")
    else:
        if not rhat_bad.empty:
            lines.append(f"✗ R̂ violations ({len(rhat_bad)} parameters):")
            lines.append(rhat_bad[["r_hat"]].to_string())
        if not ess_bulk_bad.empty:
            lines.append(f"✗ ESS-bulk violations ({len(ess_bulk_bad)} parameters):")
            lines.append(ess_bulk_bad[["ess_bulk"]].to_string())
        if not ess_tail_bad.empty:
            lines.append(f"✗ ESS-tail violations ({len(ess_tail_bad)} parameters):")
            lines.append(ess_tail_bad[["ess_tail"]].to_string())

    report = "\n".join(lines)
    print(report)

    if outpath is not None:
        full_path = os.path.join(outpath, "model_diagnostics.txt")
        with open(full_path, "w", encoding="utf-8") as f:
            f.write(report + "\n\n")
            f.write("=== Full summary ===\n")
            f.write(summary.to_string())
            f.write("\n")

    return (summary, {
        "rhat_violations": rhat_bad,
        "ess_bulk_violations": ess_bulk_bad,
        "ess_tail_violations": ess_tail_bad,
        "passed": passed,
    })


def divergence_checks(trace_negbin, outpath = None):
    divergent = trace_negbin.sample_stats["diverging"].values
    n_divergent = int(divergent.sum())
    n_total = divergent.size

    lines = [f"Divergent transitions: {n_divergent}"]
    if n_divergent == 0:
        percentage = 0
        lines.append(f"No divergent transitions (0 / {n_total} draws).")
    else:
        percentage = 100 * (n_divergent / n_total)
        lines.append(f"{n_divergent}/{n_total} ({percentage:.2f}%) divergent transitions observed.")

    report = "\n".join(lines)
    print(report)

    if outpath is not None:
        full_path = os.path.join(outpath, "model_diagnostics.txt")
        with open(full_path, "a", encoding="utf-8") as f:
                f.write("\n=== Divergence check ===\n")
                f.write(report + "\n")
        print(f"  Saved to: {full_path}")

    return (n_divergent, percentage)


def traceplots_checks(trace_negbin, outpath=None):
    figs = az.plot_trace(trace_negbin)
    figs2 = figs.ravel()[0].get_figure()
    figs2.suptitle("Trace plots — key parameters", y=1.01, fontsize=13)
    plt.tight_layout()

    if outpath is not None:
        # assumes outpath is a path
        full_path = os.path.join(outpath, "trace_plots.png")
        figs2.savefig(full_path, bbox_inches="tight", dpi=150)
        print(f"  Saved to: {full_path}")
    plt.close(figs2)
    return


# Prior Predictive checks
def prior_predictive_checks(model, inputs_dict, outpath=None):
    with model:
        prior_pred = pm.sample_prior_predictive(samples=1000, random_seed=42)

    Y = inputs_dict["Y"]
    sim_y = prior_pred.prior_predictive["Y_obs"].values.flatten()

    report = (
    "===== DATA SUMMARY STATISTICS =====\n\n"
    "== Prior Predictive Simulated Data ==\n"
    f"Min: {sim_y.min():.3f}\n"
    f"Max: {sim_y.max():.3f}\n"
    "== Real Data (Y) ==\n"
    
    f"Y min: {Y.min()}\n"
    f"Y max: {Y.max()}\n"
    "== Median Comparisons ==\n"
    f"Prior predictive median: {np.median(sim_y):.1f}\n"
    f"Real data median: {np.median(Y):.1f}\n\n"
    "    ===== INTERPRETATION =====\n"
    "Priors should be wider than the observed data - this is so to avoid ruling out plausible values.\n"
    "A prior predictive maximum around 3-5X times the observed one is informative, but loose enough not to constrain data.\n"
    "The prior predictive median should be within roughly 3× of the observed median.\n"
    "A larger gap suggests the prior is centered on the wrong scale.\n"
    )

    print(report)

    if outpath is not None:
        # assumes outpath is a path
        full_path = os.path.join(outpath, "Prior_predictive_checks.txt")
        with open(full_path, "a") as f:
            f.write(report)
        print(f"  Saved to: {full_path}")


    return


def posterior_predictive_check(model, trace_negbin, random_seed=42, outpath=None):
    with model:
        ppc = pm.sample_posterior_predictive(trace_negbin, random_seed=random_seed)

        az.plot_ppc(ppc, num_pp_samples=100)

    if outpath is not None:
        full_path = os.path.join(outpath, "posterior_predictive_check.png")
        # Save the current figure
        plt.savefig(full_path, dpi=300, bbox_inches='tight')

    plt.show()

    return


def compile_results_factor(trace_negbin, inputs_dict):
    mu_factor_logit_samples = trace_negbin.posterior['mu_factor_logit'].values  # shape: (chains, draws, n_tissues)

    # Transform to probability scale
    factor_efficiency_samples = 1 / (1 + np.exp(-mu_factor_logit_samples))  # invlogit

    # Calculate summary statistics for each tissue
    results = []

    for i, level in enumerate(inputs_dict["levels"]):
        # Get samples for this tissue (flatten chains and draws)
        samples = factor_efficiency_samples[:, :, i].flatten()

        # Calculate statistics
        mean_eff = samples.mean()
        median_eff = np.median(samples)
        ci_lower = np.percentile(samples, 2.5)
        ci_upper = np.percentile(samples, 97.5)
        ci_width = ci_upper - ci_lower

        results.append({
            'Level': level,
            'Mean RTA': mean_eff,
            'Median RTA': median_eff,
            '95% CI Lower': ci_lower,
            '95% CI Upper': ci_upper,
            'CI Width': ci_width
        })

    results = pd.DataFrame(results)

    return results


# Checks how much of the  posterior probability mass for theta_run[i] falls within the ROPE boundary
# That's the boundary where we say "it's practically zero" for the efficiency. We limit at 1%.
def flag_practical_zero_runs(trace_negbin, inputs_dict, rope_upper=0.01, mass_threshold = 0.95):
    """
    Quantify the posterior mass of theta_run within the ROPE interval [0, rope_upper]
    and provide an interpretation.

    If more than 95% of the probability masss is within the ROPE, then we are 95%
     that the efficiency is around 0.
    If instead only 5% or less, then we are more sure the value is non-zero.
    Intermediate values indicate we cannot discriminate from the data.

    trace_negbin : az.InferenceData
        Posterior samples from the fitted model.
    inputs_dict : dict
    rope_upper : float
        Upper bound of the Region of Practical Equivalence to zero.
        Default 0.01 (suitable for 'standard' and 'elevated' belief regimes).
        Consider lowering (e.g. 0.003) for the 'rare' regime ?
        mass_threshold: float (probability we want for the posterior mass to be labelled)
    """
    theta_samples = trace_negbin.posterior['theta_run'].values

    # Bins for interpretation:
    a = 1 - mass_threshold
    b = mass_threshold
    bins = [-np.inf, a, b, np.inf]
    labels = ["Efficiency credibly non-zero",
              "Inconclusive - cannot discriminate from data",
              "Efficiency too low, practically zero"
    ]
    results = []

    for i, run_name in enumerate(inputs_dict["run_names"]):
        samples = theta_samples[:, :, i].flatten()
        mass_in_rope = (samples <= rope_upper).mean()
        results.append({
            'Run': run_name,
            'P(theta <= ROPE_upper)': mass_in_rope,
        })

    results_df = pd.DataFrame(results)

    results_df['ROPE_interpretation'] = pd.cut(
        results_df['P(theta <= ROPE_upper)'],
        bins=bins,
        labels=labels,
        right=False,  # left-inclusive: [lower, upper)
    )

    return results_df


def compile_results_run(trace_negbin, inputs_dict, mass_threshold = 0.95, rope_upper = 0.01):
    run_efficiency_samples = trace_negbin.posterior['theta_run'].values  # shape: (chains, draws, n_tissues)

    # Calculate summary statistics for each tissue
    results = []

    for i, run_name in enumerate(inputs_dict["run_names"]):
        # Get samples for this tissue (flatten chains and draws)
        samples = run_efficiency_samples[:, :, i].flatten()

        # Calculate statistics
        results.append({
            'Run': run_name,
            'Mean RTA': samples.mean(),
            'Median RTA': np.median(samples),
            '95% CI Lower': np.percentile(samples, 2.5),
            '95% CI Upper': np.percentile(samples, 97.5),
            'CI Width': np.percentile(samples, 97.5) - np.percentile(samples, 2.5)
        })

    results = pd.DataFrame(results)

    df = pd.DataFrame({
        'Run': np.asarray(inputs_dict['run_names'])[inputs_dict['run_idx']],
        'level': np.asarray(inputs_dict['levels'])[inputs_dict['levels_idx']],
        'total_counts': inputs_dict["X"] + inputs_dict["Y"],
    })

    merged = pd.merge(results, df, on="Run", how="left")

    # Find null runs (where 95% prob mass for theta_run is within ROPE boundary)
    null_runs = flag_practical_zero_runs(trace_negbin = trace_negbin, inputs_dict = inputs_dict,
                                         mass_threshold=mass_threshold, rope_upper = rope_upper)
    merged = pd.merge(merged, null_runs, on="Run", how="left")

    # Flags datasets/runs where the CI width is twice the median CI width for that tissue.
    tissue_median_width = merged.groupby("level")["CI Width"].transform("median")
    merged["CI_width_ratio"] = merged["CI Width"] / tissue_median_width
    merged["Wide_CI_Flag"] = merged["CI_width_ratio"] > 2
    merged.drop(columns="CI_width_ratio", inplace=True)

    merged = merged.sort_values(by="Wide_CI_Flag", ascending=False)

    return merged


def factor_plot(df, inputs_dict, outpath=None):
    df = df.sort_values('Mean RTA', ascending=False)

    # Visualization
    fig, ax = plt.subplots(figsize=(10, max(6, inputs_dict["n_levels"] * 0.4)))

    # Sort by mean efficiency for plotting
    df_plot = df.sort_values('Mean RTA')

    y_pos = np.arange(len(df_plot))
    means = df_plot['Mean RTA'].values
    ci_lower = df_plot['95% CI Lower'].values
    ci_upper = df_plot['95% CI Upper'].values

    # Plot error bars
    ax.errorbar(means, y_pos,
                xerr=[means - ci_lower, ci_upper - means],
                fmt='o', markersize=8, capsize=5, capthick=2,
                color='steelblue', ecolor='steelblue', alpha=0.8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df_plot["Level"])
    ax.set_xlabel('RTA', fontsize=12)
    ax.set_title(f'{inputs_dict["factor"]}-Level RTA Estimates (95% Credible Intervals)', fontsize=14)
    ax.grid(True, alpha=0.3, axis='x')
    ax.axvline(x=0.5, color='gray', linestyle='--', alpha=0.5, label='50% efficiency')
    ax.legend()
    plt.tight_layout()

    if outpath is not None:
        fac = inputs_dict["factor"]
        full_path = os.path.join(outpath, f"{fac}_RTA_estimates.png")
        # Save the current figure
        plt.savefig(full_path, dpi=300, bbox_inches='tight')

    plt.show()

    return


def runs_plot(df, mass_threshold = 0.95, outpath=None):
    # Sort by level first, then by Mean RTA within each level
    df = df.sort_values(['level', 'Mean RTA']).reset_index(drop=True)

    # Assign x positions with gaps between groups
    gap = 1  # extra space between groups
    x_pos = []
    group_centers = {}
    current_x = 0
    for level, sub in df.groupby('level', sort=False):
        positions = np.arange(current_x, current_x + len(sub))
        x_pos.extend(positions)
        group_centers[level] = positions.mean()
        current_x += len(sub) + gap
    df['x_pos'] = x_pos

    # Color by level
    levels = df['level'].unique()
    cmap = plt.get_cmap('tab10')
    level_colors = {lv: cmap(i % 10) for i, lv in enumerate(levels)}

    # Assign marker shape based on ROPE_interpretation
    def resolve_marker(interp):
        if isinstance(interp, str):
            interp_lower = interp.lower()
            if 'inconclusive' in interp_lower:
                return 'x'
            if 'practically zero' in interp_lower:
                return '^'
        return 'o'  # default shape for everything else (triangle)

    df['marker'] = df['ROPE_interpretation'].apply(resolve_marker)

    # Plot
    fig, ax = plt.subplots(figsize=(max(10, len(df) * 0.35), 6))
    means = df['Mean RTA'].values
    ci_lower = df['95% CI Lower'].values
    ci_upper = df['95% CI Upper'].values

    # Plot each (level, marker) combination separately so shapes come through
    for lv in levels:
        for marker in df.loc[df['level'] == lv, 'marker'].unique():
            mask = (df['level'] == lv) & (df['marker'] == marker)
            if not mask.any():
                continue
            # Only label the "main" marker in the legend to avoid duplicates
            label = lv if marker == df.loc[df['level'] == lv, 'marker'].mode().iloc[0] else None
            ax.errorbar(df.loc[mask, 'x_pos'], df.loc[mask, 'Mean RTA'],
                        yerr=[means[mask] - ci_lower[mask], ci_upper[mask] - means[mask]],
                        fmt=marker, markersize=7, capsize=4, capthick=1.5,
                        color=level_colors[lv], ecolor=level_colors[lv],
                        alpha=0.85, label=label, linestyle='none')

    # Add shape legend entries (one per interpretation category)
    from matplotlib.lines import Line2D
    shape_handles = [
        Line2D([0], [0], marker='x', color='gray', linestyle='none',
               markersize=7, label='Inconclusive from data'),
        Line2D([0], [0], marker='^', color='gray', linestyle='none',
               markersize=7, label=f'Probably zero ({mass_threshold*100}%)'),
        Line2D([0], [0], marker='o', color='gray', linestyle='none',
               markersize=7, label=f'Likely non-zero ({mass_threshold*100}%)'),
    ]

    # Run names on x-axis
    ax.set_xticks(df['x_pos'])
    ax.set_xticklabels(df['Run'], rotation=45, ha='right', fontsize=9)

    # Level labels above the plot
    ax2 = ax.secondary_xaxis('top')
    ax2.set_xticks(list(group_centers.values()))
    ax2.set_xticklabels(list(group_centers.keys()), fontsize=11, fontweight='bold')
    ax2.tick_params(axis='x', length=0, pad=8)

    # Visual separators between groups
    boundaries = []
    current_x = 0
    for level, sub in df.groupby('level', sort=False):
        current_x += len(sub)
        boundaries.append(current_x + gap / 2 - 0.5)
        current_x += gap
    for b in boundaries[:-1]:
        ax.axvline(b, color='lightgray', linestyle=':', alpha=0.6)

    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='50% efficiency')
    ax.set_ylabel('RTA', fontsize=12)
    ax.set_title('Run-Level RTA Estimates by Group (95% Credible Intervals)', fontsize=14)
    ax.grid(True, alpha=0.3, axis='y')

    # Combined legend: level colours + shape meanings
    level_handles, level_labels = ax.get_legend_handles_labels()
    ax.legend(handles=level_handles + shape_handles,
              loc='best', fontsize=9, framealpha=0.9)

    plt.tight_layout()

    if outpath is not None:
        full_path = os.path.join(outpath, "Run_RTA_estimates.png")
        plt.savefig(full_path, dpi=300, bbox_inches='tight')

    plt.show()

    return


def run_HBM(df, metadata, coordinates, factor, outpath=None, random_state=42, min_members=10, belief = "standard", mass_threshold = 0.95, rope_upper = 0.01):
    # Keep only datasets where baseline is expressed.
    _, good_runs = test_baseline_vs_noise(df, coordinates, random_state)
    df = df[good_runs].copy()

    # Prepare input data, stores length of each slice.
    data, lengths = prepare_input(df, coordinates)

    # Merge with metadata
    data = merge_meta(data, metadata, factor)

    # NB. "scrapped" should be reported for manual inspection!
    filtered_df, scrapped = first_filter(data)

    # NB. We want at least 10 members per group, when modelling.
    filtered_2_df, scraps = second_filter(filtered_df, factor, min_members)

    # Format the data for the analysis. Stores everything in a dictionary.
    inputs_dict = format_data(filtered_2_df, factor)

    # Assemble model
    model = assemble_model(inputs_dict, lengths, belief= belief)
    trace_negbin = sample_posterior(model)

    # Run checks - Prior Predictive
    prior_predictive_checks(model, inputs_dict, outpath)
    # Run checks - Convergence
    summary, viols = check_mcmc_diagnostics(trace_negbin, outpath=outpath)
    n_divergent, percentage = divergence_checks(trace_negbin, outpath)
    # Run checks - trace plots
    traceplots_checks(trace_negbin, outpath)
    # Run checks - Posterior Predictive
    posterior_predictive_check(model, trace_negbin, random_state, outpath)

    # Compiles the table. One also flags run outliers (Wide_CI_Flag)
    res = compile_results_factor(trace_negbin, inputs_dict)
    runs = compile_results_run(trace_negbin = trace_negbin, inputs_dict = inputs_dict, mass_threshold = mass_threshold, rope_upper = rope_upper)

    # Save res and runs
    full_path = os.path.join(outpath, f"{factor}_RTA_estimates.csv")
    res.to_csv(full_path, index=True)

    full_path = os.path.join(outpath, f"Run_RTA_estimates.csv")
    runs.to_csv(full_path, index=True)

    # Makes plots and saves
    factor_plot(res, inputs_dict, outpath)
    runs_plot(runs, mass_threshold = mass_threshold, outpath = outpath)

    return res, runs