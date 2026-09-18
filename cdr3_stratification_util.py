"""Helpers for testing whether selection at IGHV sites depends on the CDR-H3 context.

The question is whether entrenchment at sites 33, 35, 50 and 52 is genuinely driven by
contact with the CDR-H3 loop. One way to find out is to stratify clonal families by CDR-H3
length, charge or hydrophobicity. CDR-H3 is the one genetically
uncoupled partner inside DASM's input sequence, so unlike the light chain or the antigen it
can be varied while holding the germline background fixed.

One design choice shapes what is here: the tests these helpers feed are run
*relatively*, comparing CDR-H3 sensitivity at contacting against non-contacting sites,
because the absolute version is confounded by donor, J gene and SHM load, all of which
co-vary with CDR-H3 composition and would shift every site alike.

Effective sample size is the number of clonal families, not parent-child pairs: all PCPs in
a family share one CDR-H3, so anything built on these helpers should treat the family as
the unit.

This module holds only shared data preparation. The statistics live in the notebooks:
cdr3_dependence_of_selection.ipynb and cdr3_pairing_bias_with_entrenched_sites.ipynb.
"""

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# The CDR-H3 descriptors we stratify on.
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
    supplies the counterexample -- are whatever the data says they are.

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
