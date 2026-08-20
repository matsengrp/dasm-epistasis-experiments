"""Helpers for testing whether selection at IGHV sites depends on the CDR-H3 context.

Reviewer 3 of the MBE submission asks whether entrenchment at sites 33, 35, 50 and 52 is
genuinely driven by contact with the CDR-H3 loop, and proposes stratifying clonal families
by CDR-H3 length, charge or hydrophobicity to find out. CDR-H3 is the one genetically
uncoupled partner inside DASM's input sequence, so unlike the light chain or the antigen it
can be varied while holding the germline background fixed.

Two design choices shape what is here. First, we run the test *relatively* -- comparing
CDR-H3 sensitivity at contacting versus non-contacting sites -- because the absolute
version is confounded by donor, J gene and SHM load, all of which co-vary with CDR-H3
composition and would shift every site alike. Second, every DASM-based test has a
model-free twin built from observed substitution counts, so the answer does not rest on
trusting DASM.

Effective sample size is the number of clonal families, not parent-child pairs: all PCPs in
a family share one CDR-H3. Standard errors are clustered on family throughout, which
matters a lot -- on a cluster-constant predictor like these, ignoring it inflates
significance about four-fold.
"""

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# statsmodels.api and statsmodels.formula.api are broken in netam_env (statsmodels 0.14.4
# against scipy 1.16 raises ImportError on scipy._lib._util._lazywhere). The estimator
# classes themselves import fine, so we use those directly and build design matrices with
# pd.get_dummies instead of going through the formula interface.
from statsmodels.genmod.families import Poisson
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.regression.linear_model import OLS, WLS

# The descriptors Reviewer 3 names.
CDR3_DESCRIPTORS = ["cdr3_length", "cdr3_charge", "cdr3_hydrophobicity"]

# The V gene encodes up to site 94 in Chothia numbering; 95-102 is the CDR-H3 itself and
# 103 onwards is J-encoded FR4. Restricting to sites <= 94 keeps the analysis inside the
# region entrenchment is defined on, and incidentally avoids sites whose apparent CDR-H3
# burial is chain connectivity to the deleted window rather than a long-range contact.
MAX_V_REGION_SITE = 94


def site_number(site):
    """Numeric part of a Chothia site label, so insertion codes like 52A sort correctly."""
    digits = "".join(c for c in str(site) if c.isdigit())
    return int(digits) if digits else -1


def v_region_only(df, site_col="site"):
    """Restrict to V-gene-encoded sites (<= 94)."""
    sites = df.index if site_col is None else df[site_col]
    return df[[site_number(s) <= MAX_V_REGION_SITE for s in sites]]


# A site counts as CDR-H3-contacting when the CDR-H3 reaches it in most structures.
# Contact frequency is used rather than median burial because the per-structure burial
# distributions are bimodal -- site 52 has a median of 0.004 yet is contacted in half of
# structures, and site 2 is contacted in 98% but shallowly. Frequency is also the quantity
# that matches what a DASM selection factor averages over: a repertoire of many CDR-H3s.
CDR3_CONTACT_THRESHOLD = 0.5


def add_cdr3_descriptors(pcp_df):
    """Add per-PCP CDR-H3 sequence and physicochemical descriptors.

    The CDR-H3 span comes from cdr3_codon_start/cdr3_codon_end on the *parent* sequence,
    which excludes the conserved flanking Cys and Trp. Reading it off the parent rather
    than the child makes the descriptor a property of the branch's starting point, matching
    how the selection factor is defined.

    A span is valid when it is annotated, runs forwards, and is a whole number of codons.
    Pairs failing that get a NaN CDR-H3 and drop out of every fit downstream. On
    v1rodriguez exactly 4 of 21752 pairs fail, all of them the one clonal family whose
    annotation places the CDR-H3 end three nucleotides before its start; the count is
    returned to the caller so a dataset where this is not negligible cannot slip by.

    Nothing guards against stop or ambiguous codons in the CDR-H3. These are productive
    repertoires and contain none (verified: 0 of 21752), so if one ever appeared the data
    would not be what we think it is, and letting ProteinAnalysis raise is the right
    outcome.
    """
    df = pcp_df.copy()
    start = pd.to_numeric(df.cdr3_codon_start, errors="coerce")
    end = pd.to_numeric(df.cdr3_codon_end, errors="coerce")
    valid = start.notna() & end.notna() & (end > start) & ((end - start) % 3 == 0)

    cdr3_aa = [
        str(Seq(parent[int(s) : int(e)]).translate()) if ok else None
        for ok, s, e, parent in zip(valid, start, end, df.parent)
    ]

    df["cdr3_aa"] = cdr3_aa
    seqs = df.cdr3_aa.dropna()
    analyses = seqs.apply(ProteinAnalysis)

    df["cdr3_length"] = seqs.str.len()
    # Biopython's charge_at_pH solves Henderson-Hasselbalch over published side-chain pKa
    # values, rather than the usual shortcut of counting K/R minus D/E. It also adds free
    # N- and C-terminal charges, which the CDR-H3 does not have inside an intact chain;
    # that contribution is near-constant across sequences (CDR-H3s begin with the same
    # residues after the conserved Cys) and the descriptor is standardised before fitting,
    # so it shifts the intercept rather than the slope.
    df["cdr3_charge"] = analyses.apply(lambda a: a.charge_at_pH(7.0))
    # GRAVY: mean Kyte-Doolittle hydropathy (Kyte & Doolittle 1982, J Mol Biol 157:105).
    df["cdr3_hydrophobicity"] = analyses.apply(lambda a: a.gravy())
    return df


def classify_sites(burial_df, numbering_scheme="chothia"):
    """Label every site by entrenchment status and CDR-H3 contact, both from data.

    The two axes are orthogonal and were previously conflated in a hand-written panel,
    which mislabelled site 94 as non-entrenched when it is entrenched between families.

    Entrenchment comes from the saved entrenchment calls, contact from structural burial.
    Nothing here is hand-picked, so the non-entrenched contacting sites -- the class that
    supplies Reviewer 3's counterexample -- are whatever the data says they are.

    Note that a site can be entrenched at both levels; ``entrenchment`` reports the
    within-family call where both apply, since that is the paper's primary analysis.
    """
    from utils import load_entrenched_sites

    _, _, _, _, within, between = load_entrenched_sites(numbering_scheme)
    within_sites = set(within.site.astype(str))
    between_sites = set(between.site.astype(str))

    out = burial_df.copy()
    out["entrenchment"] = [
        "within-family" if s in within_sites
        else "between-family" if s in between_sites
        else "not entrenched"
        for s in out.index
    ]
    out["cdr3_contact"] = out.cdr3_contact_freq > CDR3_CONTACT_THRESHOLD
    out["site_class"] = out.entrenchment + np.where(
        out.cdr3_contact, ", CDR-H3 contact", ", no contact"
    )
    return out


def cdr3_burial_by_site(
    sasa_path="_output/sasa_human_chothia_anarci.csv",
    v_families=("IGHV1", "IGHV3"),
    min_structures=30,
):
    """Change in RSA per site when the CDR-H3 is deleted, over SAbDab structures.

    Uses ``rsa_vdj_junction_effect``, the quantity the paper plots and analyses: its CDR3
    window is positions 95-102, one position beyond the Chothia CDR-H3 boundaries on each
    side, capturing the variable part of the V(D)J junction. The similarly named
    ``rsa_cdr3_effect`` column applies the strict Chothia boundaries and is *not* what
    solvent_accessibility_analysis.ipynb uses. The two differ substantially -- site 2 is
    0.009 under one and 0.171 under the other -- so the name is not a safe guide.

    Sign flipped so larger means more buried. ``cdr3_contact_freq`` is the fraction of
    structures with any contact, which is a distinct quantity from the median depth: a site
    can be contacted rarely but deeply.
    """
    cols = [
        "pdb_id", "v_gene_heavy", "site", "chain_id", "heavy_chain_id",
        "rsa_vdj_junction_effect", "rsa_antigen_effect", "rsa_light_effect",
    ]
    d = pd.read_csv(sasa_path, dtype={"site": str}, usecols=cols)
    d = d[d.chain_id == d.heavy_chain_id]
    d = d[d.v_gene_heavy.str.extract(r"(IGHV\d+)")[0].isin(v_families)]
    d = d.dropna(subset=["rsa_vdj_junction_effect"])

    g = d.groupby("site").agg(
        cdr3_burial=("rsa_vdj_junction_effect", lambda x: -x.median()),
        cdr3_contact_freq=("rsa_vdj_junction_effect", lambda x: float((x < -1e-9).mean())),
        antigen_burial=("rsa_antigen_effect", lambda x: -x.mean()),
        light_burial=("rsa_light_effect", lambda x: -x.mean()),
        n_structures=("pdb_id", "nunique"),
    )
    return g[g.n_structures >= min_structures]


def site_level_lsf(dasm_df):
    """One selection value per (PCP, site): the mean log selection factor over 1nt targets.

    Measures how disfavoured *any* departure from the germline residue is, rather than one
    specific reciprocal substitution. Much better powered than the pair-level readout.
    """
    return (
        dasm_df.groupby(["pcp_index", "site", "v_gene", "parent_aa"], observed=True)
        .log_selection_factor.mean()
        .reset_index(name="mean_lsf")
    )


def _design(d, factor_cols):
    """Intercept, standardised descriptor 'x', and dummy codes for each factor."""
    parts = [pd.Series(1.0, index=d.index, name="const"), d.x]
    for col in factor_cols:
        if d[col].nunique() > 1:
            parts.append(
                pd.get_dummies(d[col].astype(str), prefix=col, drop_first=True).astype(float)
            )
    return pd.concat(parts, axis=1)


def _standardise(d, descriptor):
    """Centre and scale the descriptor so slopes read as 'per one standard deviation'."""
    d = d.copy()
    d["x"] = (d[descriptor] - d[descriptor].mean()) / d[descriptor].std()
    return d
def add_branch_exposure(site_sub_df, pseudocount=0.5):
    """Per-branch neutral exposure: synonymous mutations per site, as a Poisson offset.

    Synonymous divergence is the clock the paper's rate analysis already uses, since it
    measures elapsed mutation without being shaped by amino-acid selection.

    The pseudocount matters. About one branch in six carries no synonymous mutation, and
    log(0) is undefined, so a bare synonymous frequency would force those branches out of
    the model -- biasing toward heavily mutated branches even though a short branch can
    still carry a substitution. Adding 0.5 is the same smoothing rates_analysis_util
    applies to its counts.
    """
    synonymous = np.where(site_sub_df.parent_aa == site_sub_df.child_aa,
                          site_sub_df.nucleotide_mutation_count, 0)
    per_branch = pd.Series(synonymous, index=site_sub_df.index).groupby(
        site_sub_df.pcp_index).transform("sum")
    n_codons = site_sub_df.groupby("pcp_index").pcp_index.transform("size")
    return site_sub_df.assign(branch_exposure=(per_branch + pseudocount) / (n_codons * 3))


def observed_substitution_counts(site_sub_df, remove_leaves=True):
    """Per-PCP indicator of a non-synonymous substitution: the model-free readout.

    Keeps only branches whose parent still carries the germline codon, so a substitution is
    a genuine departure from the germline state, and drops leaf branches following the
    paper's rate analysis.
    """
    d = site_sub_df
    if remove_leaves:
        d = d[~d.child_is_leaf]
    d = d[d.is_germline_codon == True]  # noqa: E712
    return d.assign(
        substituted=((d.nucleotide_mutation_count > 0) & (d.child_aa != d.parent_aa)).astype(int)
    )


def fit_cdr3_rate_slope(df, descriptor, strata=("v_gene",), min_events=15):
    """Poisson regression of observed substitutions on a CDR-H3 descriptor. No DASM.

    A continuous fit is used rather than Reviewer 3's tercile comparison because entrenched
    sites are by construction the sites with fewest substitutions -- splitting them three
    ways leaves single-digit counts. Terciles are still produced for plotting.

    Dummy levels with no substitution events are dropped first. They are perfectly
    separated, so their coefficients run to minus infinity and IRLS fails to converge, even
    though the coefficient of interest is well identified.
    """
    d = df.dropna(subset=[descriptor, "branch_exposure", "family_key"])
    if d.substituted.sum() < min_events or d[descriptor].std() == 0:
        return None
    d = _standardise(d, descriptor)
    X = _design(d, list(strata))
    separated = [c for c in X.columns if c not in ("const", "x")
                 and d.substituted[X[c] > 0].sum() == 0]
    X = X.drop(columns=separated)

    fit = GLM(d.substituted, X, family=Poisson(), offset=np.log(d.branch_exposure)).fit(
        cov_type="cluster", cov_kwds={"groups": d.family_key}
    )
    return {
        "descriptor": descriptor,
        "log_rate_ratio": fit.params["x"],
        "se": fit.bse["x"],
        "pvalue": fit.pvalues["x"],
        "n_pcp": len(d),
        "n_events": int(d.substituted.sum()),
        "n_families": d.family_key.nunique(),
    }


def slopes_by_site(df, descriptor, fitter, **kwargs):
    """Run ``fitter`` independently at every site. Returns a tidy DataFrame.

    Sites where the fit is degenerate (singular design, non-convergence) are skipped rather
    than aborting the sweep.
    """
    rows, failed = [], []
    for site, sub in df.groupby("site", observed=True):
        try:
            res = fitter(sub, descriptor, **kwargs)
        except (np.linalg.LinAlgError, PerfectSeparationError) as e:
            failed.append(f"{site} ({type(e).__name__})")
            continue
        if res is not None:
            rows.append({"site": site, **res})
    if failed:
        print(f"  {descriptor}: no fit at {len(failed)} sites: {', '.join(failed)}")
    return pd.DataFrame(rows)
def germline_selection_strength(numbering_scheme="chothia", v_families=("IGHV1", "IGHV3"),
                                min_directions=4):
    """Per-site strength and direction of selection on germline<->germline substitutions.

    Entrenchment needs two things: at least two germline amino acids at the site to compare
    reciprocally, and purifying selection on the substitutions between them. Sites can make
    CDR-H3 contact and still fail the second condition, which is what this quantifies.

    Reads the reciprocal comparison tables written by v_families_entrenchment_dasm.ipynb.
    Each pair appears twice there, once per orientation, so log_selection_factor_1 across
    all rows covers every direction exactly once.

    Returns median log selection factor and the fraction of directions that are *favoured*
    (selection factor above the neutral rate) -- the signature of diversifying rather than
    purifying selection.
    """
    rows = []
    for family in v_families:
        d = pd.read_csv(
            f"_output/entrenchment_analysis/{numbering_scheme}/comparison_within_{family}.csv",
            dtype={"site": str},
        )
        for site, g in d.groupby("site"):
            if len(g) < min_directions:
                continue
            lsf = g.log_selection_factor_1
            rows.append({
                "v_family": family,
                "site": site,
                "n_germline_aa": g.parent_aa_1_and_target_aa_2.nunique(),
                "n_directions": len(lsf),
                "median_lsf": lsf.median(),
                "frac_favored": (lsf > 0).mean(),
            })
    return pd.DataFrame(rows)


def family_level_selection(site_lsf):
    """One row per (site, V gene, clonal family), carrying its parent-child pair count.

    Binning has to happen at the family level, because a family has one CDR-H3 and so
    cannot be split across two CDR-H3 terciles. Collapsing costs nothing here: at depth 2 a
    family contributes one or two pairs, and where there are two they are two children of
    the *same* parent node with an identical parent sequence. A DASM selection factor is a
    function of the parent sequence alone, so those pairs carry identical selection factors
    -- verified on 399 of 399 multi-pair families checked, every value equal. The collapse
    is therefore exact, not an averaging choice, and ``n_pcp`` is retained so the pair
    counts stay comparable with the rest of the analysis.

    This does *not* carry over to the observed-substitution analysis in section 8. There the
    two pairs have different children and so different substitutions, and both are kept.
    """
    keys = ["site", "v_family", "v_gene", "parent_aa", "family_key"]
    agg = {d: (d, "first") for d in CDR3_DESCRIPTORS}
    return site_lsf.groupby(keys, observed=True).agg(
        value=("mean_lsf", "first"), n_pcp=("mean_lsf", "size"), **agg
    ).reset_index()
def entrenched_units(numbering_scheme="chothia"):
    """(V family, site, germline residue) triples that participate in an entrenched pair."""
    from utils import load_entrenched_sites

    _, aas, _, _, _, _ = load_entrenched_sites(numbering_scheme)
    return {(r.v_family, str(r.site), r.amino_acid) for r in aas.itertuples(index=False)}
def partial_r2_by_site(family_df, descriptor, min_families=60, n_bins=None):
    """Fraction of within-V-gene selection variance a CDR-H3 property explains, per site.

    This is the statistic the section 6 test reports, and it is the one that answers the
    question as posed: holding V gene -- and therefore germline residue -- fixed, how much
    does the CDR-H3 property predict the selection factor here?

    It is reported instead of a regression slope or a median shift because both of those
    carry the units of the outcome, and the spread of log selection factors differs from
    site to site (interquartile range is 1.26x larger at CDR-H3-contacting sites). A slope
    of a given size therefore means something different at different sites, and comparing
    slopes across sites measures the spread as much as the dependence. Partial R-squared is
    a ratio of variances, so the site's scale cancels.

    It also folds in exactly what a stratified median comparison leaves out: the numerator
    responds to separation between CDR-H3 groups, the denominator to spread within them.

    ``n_bins`` controls the functional form. The default fits the descriptor as a single
    linear term, which assumes selection responds monotonically. Setting ``n_bins`` instead
    splits it into quantile bins and fits dummies, which captures any shape -- a threshold,
    a U, a plateau -- at the cost of n_bins-1 parameters instead of one. That extra cost
    inflates R-squared mechanically, so binned and linear values are not comparable to each
    other; only binned-versus-binned across sites is.

    Computed at the family level, since a family has one CDR-H3.
    """
    rows = []
    for site, g in family_df.groupby("site", observed=True):
        g = g.dropna(subset=[descriptor, "value"])
        if len(g) < min_families or g[descriptor].std() == 0:
            continue
        dummies = pd.get_dummies(g.v_gene.astype(str), drop_first=True).astype(float)
        base = pd.concat([pd.Series(1.0, index=g.index, name="const"), dummies], axis=1)
        if n_bins is None:
            extra = g[descriptor].rename("x").to_frame()
        else:
            binned = pd.qcut(g[descriptor], n_bins, labels=False, duplicates="drop")
            if pd.Series(binned).nunique() < 2:
                continue
            extra = pd.get_dummies(pd.Series(binned, index=g.index),
                                   prefix="bin", drop_first=True).astype(float)
        full = pd.concat([base, extra], axis=1)
        ssr_base = OLS(g.value, base).fit().ssr
        ssr_full = OLS(g.value, full).fit().ssr
        rows.append({
            "site": site,
            "descriptor": descriptor,
            "partial_r2": (ssr_base - ssr_full) / ssr_base,
            "n_terms": extra.shape[1],
            "n_families": len(g),
            "n_pcp": int(g.n_pcp.sum()),
            "n_vgenes": g.v_gene.nunique(),
        })
    return pd.DataFrame(rows)


def r2_by_unit(family_df, descriptor, min_families=60, n_bins=None,
               group_cols=("site", "v_family", "v_gene", "parent_aa"),
               numbering_scheme="chothia"):
    """R-squared of a CDR-H3 property within each (site, V gene) unit.

    The per-site version pools all V genes at a site and fits one shared slope, which
    dilutes any effect specific to particular germline configurations: at a single site some
    V genes carry an entrenched germline residue and others do not. Here each V gene at each
    site is its own unit, which is the condition Reviewer 3 actually states -- same V gene,
    same germline residue -- and lets entrenched and non-entrenched units be compared
    directly.

    The germline residue needs no separate grouping. Rows are filtered to parent codon equal
    to germline codon, so within a (site, V gene) the germline residue is fixed by
    construction.

    With no covariates left inside a unit, this is the ordinary R-squared of selection on
    the descriptor. With ``n_bins`` unset that is the squared Pearson correlation, which
    detects only a monotone linear relationship. Setting ``n_bins`` fits bin dummies
    instead, giving the correlation ratio eta-squared: the one-way ANOVA effect size, which
    responds to any difference between bins regardless of shape.

    Note eta-squared is not guaranteed to exceed r-squared. A step function cannot
    reproduce a straight line, so coarse bins can capture *less* of a genuinely linear
    relationship than a single linear term does.
    """
    entrenched = entrenched_units(numbering_scheme)
    rows = []
    group_cols = list(group_cols)
    for key, g in family_df.groupby(group_cols, observed=True):
        keyd = dict(zip(group_cols, key if isinstance(key, tuple) else (key,)))
        site = keyd.get("site")
        v_family = keyd.get("v_family", "")
        v_gene = keyd.get("v_gene", "(pooled)")
        parent_aa = keyd.get("parent_aa", "(pooled)")
        g = g.dropna(subset=[descriptor, "value"])
        if len(g) < min_families or g[descriptor].std() == 0 or g.value.std() == 0:
            continue
        if n_bins is None:
            r = np.corrcoef(g[descriptor], g.value)[0, 1]
            r2, n_terms = r**2, 1
        else:
            binned = pd.qcut(g[descriptor], n_bins, labels=False, duplicates="drop")
            if pd.Series(binned).nunique() < 2:
                continue
            grand = g.value.mean()
            ss_total = ((g.value - grand) ** 2).sum()
            ss_between = sum(
                (binned == b).sum() * (g.value[binned == b].mean() - grand) ** 2
                for b in np.unique(binned)
            )
            r, r2 = np.nan, ss_between / ss_total
            n_terms = int(pd.Series(binned).nunique()) - 1
        rows.append({
            "site": site, "v_family": v_family, "v_gene": v_gene,
            "germline_aa": parent_aa, "descriptor": descriptor,
            "r2": r2, "r": r, "n_terms": n_terms, "n_families": len(g),
            "entrenched": any(
                (vf, site, aa) in entrenched
                for vf in ([v_family] if v_family else g.v_family.unique())
                for aa in ([parent_aa] if parent_aa != "(pooled)"
                           else g.parent_aa.unique())
            ),
        })
    return pd.DataFrame(rows)


def substitution_level_selection(dasm_df):
    """One row per (family, site, V gene, germline residue, target residue).

    Keeps each substitution separate instead of averaging over the targets reachable from
    the germline codon. That average is not comparable across units: the codon determines
    which residues are one nucleotide away, so two units with the same germline residue but
    different codons average over different target sets.

    Pairs within a family are deduplicated rather than averaged -- their parent sequence is
    identical, so their selection factors are too.
    """
    keys = ["site", "v_family", "v_gene", "parent_aa", "selection_factor_target_aa",
            "family_key"]
    return (dasm_df.sort_values("pcp_index").drop_duplicates(keys)
            [keys + ["log_selection_factor"] + CDR3_DESCRIPTORS])


SUBSTITUTION_GROUPINGS = {
    # V gene matched: the germline codon is fixed, so every family in a unit is scored on
    # the identical substitution from the identical codon.
    "site + V gene": (["site", "v_family", "v_gene", "parent_aa",
                       "selection_factor_target_aa"], ["site", "v_gene"]),
    # Pooled over V genes sharing a germline residue. All families are still scored on the
    # same residue-to-residue substitution, but the codon and the V-gene background may
    # differ, and only V genes whose codon can reach the target contribute.
    "site + germline residue": (["site", "parent_aa", "selection_factor_target_aa"],
                                ["site", "parent_aa"]),
}


def r2_by_substitution(sub_df, descriptor, n_bins=3, min_families=60,
                       grouping="site + V gene", numbering_scheme="chothia"):
    """Variance in selection explained by a CDR-H3 property, per specific substitution.

    Each unit is one (site, V gene, germline residue -> target residue). Because the
    germline codon is fixed by the V gene under the germline-codon filter, every family in
    a unit is scored on the *same* substitution, which is what makes units comparable.

    With ``n_bins`` set, the CDR-H3 property is binned within the unit's family pool and the
    statistic is the correlation ratio eta-squared: between-bin sum of squares over total,
    i.e. a one-way ANOVA effect size, which responds to any shape of relationship. With
    ``n_bins=None`` it is the squared Pearson correlation instead, which detects only a
    monotone linear response and is reported as a sensitivity check.

    Entrenchment is defined per amino-acid pair, so units are labelled against the
    (V family, site, residue, target) calls directly rather than by site.
    """
    from utils import load_entrenched_sites

    _, aas, _, _, _, _ = load_entrenched_sites(numbering_scheme)
    entrenched_pairs = {
        (r.v_family, str(r.site), r.amino_acid, r.target_amino_acid)
        for r in aas.itertuples(index=False)
    }

    unit, bin_within = SUBSTITUTION_GROUPINGS[grouping]
    unit, bin_within = list(unit), list(bin_within)
    d = sub_df.dropna(subset=[descriptor, "log_selection_factor"]).copy()

    if n_bins is None:
        # squared Pearson correlation per unit, from grouped sums
        d["_xy"] = d[descriptor] * d.log_selection_factor
        d["_xx"] = d[descriptor] ** 2
        d["_yy"] = d.log_selection_factor ** 2
        g = d.groupby(unit, observed=True)
        agg = g.agg(n_families=(descriptor, "size"), sx=(descriptor, "sum"),
                    sy=("log_selection_factor", "sum"), sxy=("_xy", "sum"),
                    sxx=("_xx", "sum"), syy=("_yy", "sum"))
        num = (agg.n_families * agg.sxy - agg.sx * agg.sy) ** 2
        den = ((agg.n_families * agg.sxx - agg.sx ** 2)
               * (agg.n_families * agg.syy - agg.sy ** 2))
        out = agg[["n_families"]].copy()
        out["r2"] = (num / den).where(den > 0)
    else:
        # bin families once per pool: the property is a family attribute, so every
        # substitution from that pool shares the same bin assignment
        fam = d[bin_within + ["family_key", descriptor]].drop_duplicates()
        fam["bin"] = fam.groupby(bin_within, observed=True)[descriptor].transform(
            lambda x: pd.qcut(x, n_bins, labels=False, duplicates="drop"))
        d = d.merge(fam[bin_within + ["family_key", "bin"]],
                    on=bin_within + ["family_key"], how="left")
        d = d.dropna(subset=["bin"])

        tot = d.groupby(unit, observed=True).log_selection_factor.agg(
            n_families="size", grand="mean", var=lambda x: x.var(ddof=0))
        per_bin = d.groupby(unit + ["bin"], observed=True).log_selection_factor.agg(
            ["size", "mean"])
        joined = per_bin.join(tot[["grand"]], on=unit)
        ss_between = ((joined["mean"] - joined["grand"]) ** 2 * joined["size"]).groupby(
            level=list(range(len(unit)))).sum()

        out = tot.copy()
        out["ss_total"] = out["var"] * out["n_families"]
        out["r2"] = (ss_between / out.ss_total).where(out.ss_total > 0)

    out = out[(out.n_families >= min_families) & out.r2.notna()].reset_index()
    out["substitution"] = out.parent_aa + "->" + out.selection_factor_target_aa
    if "v_family" in out.columns:
        out["entrenched"] = [
            (vf, s, a, t) in entrenched_pairs
            for vf, s, a, t in zip(out.v_family, out.site, out.parent_aa,
                                   out.selection_factor_target_aa)
        ]
    else:
        # pooled over V families: entrenched if the pair is called in any of them
        by_pair = {(s, a, t) for _, s, a, t in entrenched_pairs}
        out["entrenched"] = [
            (s, a, t) in by_pair
            for s, a, t in zip(out.site, out.parent_aa,
                               out.selection_factor_target_aa)
        ]
    out["descriptor"] = descriptor
    out["grouping"] = grouping
    keep = [c for c in ["site", "v_family", "v_gene", "parent_aa",
                        "selection_factor_target_aa", "substitution", "descriptor",
                        "grouping", "r2", "n_families", "entrenched"] if c in out.columns]
    return out[keep]
