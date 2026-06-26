"""
Germline misassignment risk analysis at entrenched sites.

The concern (Dunn-Walters): if two V gene alleles are very close in sequence,
an early SHM mutation at one of their differing sites could flip the V gene
assignment. If that differing site is an entrenched site, we'd mislabel a
mutated site as "germline/unmutated."

Definition of "confusable": two alleles with total nucleotide distance <= 2.
After a single mutation at the entrenched site, the sequence is equidistant
or closer to the wrong allele, so partis could pick either one.

Analysis:
1. Find all confusable V gene pairs (total nt distance <= 2)
2. For each pair, identify which sites differ and whether any are entrenched
3. Per entrenched site and V family, compute what % of V genes are "at risk"
   (have at least one confusable partner differing at that site)
4. In the Rodriguez dataset, compute what % of sequences use at-risk V genes

Usage:
    python confusable_allele_analysis.py
"""

from itertools import combinations

import pandas as pd

GERMLINE_CODONS_PATH = "germline/germline_codons_chothia.csv"
ENTRENCHED_DIR = "_output/entrenchment_analysis/chothia"
PCP_PATH = ("/home/nharel/re/dnsm-experiments-1/dasm-train/_ignore/test_output/"
            "dasm_4m-v1jaffeCC+v1tangCC-joint-ON-v1rodriguez-chothia-pcp_df.csv")

WITHIN_FAMILIES = ["IGHV1", "IGHV3", "IGHV4"]

# Two alleles are "confusable" if their total nt distance is at most this.
# At distance <= 2, a single mutation at a differing site makes the sequence
# equidistant or closer to the wrong allele.
CONFUSABLE_NT_THRESHOLD = 2


def load_germline_codons():
    return pd.read_csv(GERMLINE_CODONS_PATH, dtype={"site": str})


def load_within_family_entrenched_sites():
    """Load within-family entrenched sites: {v_family: set of site strings}."""
    entrenched_sites = {}
    for fam in WITHIN_FAMILIES:
        path = f"{ENTRENCHED_DIR}/entrenched_aa_sites_within_{fam}.csv"
        df = pd.read_csv(path, dtype={"site": str})
        entrenched_sites[fam] = set(df["site"].unique())
    return entrenched_sites


def pairwise_site_differences(codon_wide, gene1, gene2):
    """Return {site: (codon1, codon2)} for sites where the two genes differ."""
    row1 = codon_wide.loc[gene1]
    row2 = codon_wide.loc[gene2]
    diffs = {}
    for site in codon_wide.columns:
        c1, c2 = row1[site], row2[site]
        if pd.isna(c1) or pd.isna(c2):
            continue
        if c1 != c2:
            diffs[site] = (c1, c2)
    return diffs


def codon_nt_distance(c1, c2):
    return sum(a != b for a, b in zip(c1, c2))


def main():
    codons_df = load_germline_codons()
    entrenched_sites_by_family = load_within_family_entrenched_sites()

    codon_wide = codons_df.pivot(index="v_gene", columns="site", values="codon")
    aa_wide = codons_df.pivot(index="v_gene", columns="site", values="amino_acid")
    gene_to_family = (
        codons_df[["v_gene", "v_family"]]
        .drop_duplicates()
        .set_index("v_gene")["v_family"]
    )

    # ── Step 1: find all confusable pairs ──
    all_genes = sorted(codon_wide.index)
    confusable_pairs = []

    for gene1, gene2 in combinations(all_genes, 2):
        site_diffs = pairwise_site_differences(codon_wide, gene1, gene2)
        total_nt_diff = sum(codon_nt_distance(c1, c2) for c1, c2 in site_diffs.values())

        if total_nt_diff > CONFUSABLE_NT_THRESHOLD:
            continue

        fam1 = gene_to_family.get(gene1, "")
        fam2 = gene_to_family.get(gene2, "")

        confusable_pairs.append({
            "gene1": gene1,
            "gene2": gene2,
            "family1": fam1,
            "family2": fam2,
            "total_nt_distance": total_nt_diff,
            "differing_sites": site_diffs,
        })

    print("=" * 80)
    print(f"CONFUSABLE V GENE PAIRS (total nt distance <= {CONFUSABLE_NT_THRESHOLD})")
    print("=" * 80)
    print(f"\nFound {len(confusable_pairs)} confusable pairs\n")

    # ── Step 2: for each pair, which sites differ, and are any entrenched? ──
    # Build: for each (site, gene), is this gene "at risk" at this site?
    # A gene is at risk at a site if it has a confusable partner that differs at that site.
    at_risk_genes_per_site = {}  # (entrenched_site, v_family) -> set of gene names

    print("CONFUSABLE PAIRS AND THEIR DIFFERING SITES:")
    print(f"{'─' * 80}")

    for pair in confusable_pairs:
        gene1, gene2 = pair["gene1"], pair["gene2"]
        site_diffs = pair["differing_sites"]
        differing_sites = set(site_diffs.keys())

        # Check against entrenched sites for both genes' families
        fam1, fam2 = pair["family1"], pair["family2"]
        relevant_entrenched = set()
        for fam in [fam1, fam2]:
            if fam in entrenched_sites_by_family:
                relevant_entrenched |= entrenched_sites_by_family[fam]

        entrenched_diffs = sorted(differing_sites & relevant_entrenched)
        non_entrenched_diffs = sorted(differing_sites - relevant_entrenched)

        # Build detail strings
        details = []
        for site in sorted(differing_sites):
            c1, c2 = site_diffs[site]
            aa1 = aa_wide.loc[gene1, site] if site in aa_wide.columns else "?"
            aa2 = aa_wide.loc[gene2, site] if site in aa_wide.columns else "?"
            nt_d = codon_nt_distance(c1, c2)
            is_ent = site in relevant_entrenched
            marker = " **ENTRENCHED**" if is_ent else ""
            aa_str = f"{aa1}->{aa2}" if aa1 != aa2 else f"{aa1}(syn)"
            details.append(f"    site {site}: {aa_str} ({nt_d}nt, {c1}->{c2}){marker}")

        print(f"\n  {gene1}  vs  {gene2}  (dist={pair['total_nt_distance']})")
        for d in details:
            print(d)

        # Record at-risk genes
        for site in entrenched_diffs:
            for fam in [fam1, fam2]:
                if fam in entrenched_sites_by_family and site in entrenched_sites_by_family[fam]:
                    key = (site, fam)
                    if key not in at_risk_genes_per_site:
                        at_risk_genes_per_site[key] = set()
                    # Both genes in the pair are at risk
                    if gene_to_family.get(gene1) == fam:
                        at_risk_genes_per_site[key].add(gene1)
                    if gene_to_family.get(gene2) == fam:
                        at_risk_genes_per_site[key].add(gene2)

    # ── Step 3: per entrenched site and V family, % of genes at risk ──
    print(f"\n\n{'=' * 80}")
    print("AT-RISK V GENES PER ENTRENCHED SITE")
    print(f"{'=' * 80}")
    print(f"(genes that have a confusable partner differing at that entrenched site)\n")

    for v_family in WITHIN_FAMILIES:
        family_genes = set(gene_to_family[gene_to_family == v_family].index)
        n_family_genes = len(family_genes)
        ent_sites = sorted(entrenched_sites_by_family.get(v_family, set()))

        print(f"\n{v_family} ({n_family_genes} alleles in OGRDB):")

        for site in ent_sites:
            key = (site, v_family)
            risk_genes = at_risk_genes_per_site.get(key, set())
            n_risk = len(risk_genes)
            pct = 100 * n_risk / n_family_genes if n_family_genes > 0 else 0

            if risk_genes:
                gene_list = ", ".join(sorted(risk_genes))
                print(f"  Site {site}: {n_risk}/{n_family_genes} alleles at risk ({pct:.1f}%)")
                print(f"    {gene_list}")
            else:
                print(f"  Site {site}: 0/{n_family_genes} alleles at risk (0%)")

    # ── Step 4: Rodriguez dataset exposure ──
    print(f"\n\n{'=' * 80}")
    print("RODRIGUEZ DATASET EXPOSURE")
    print(f"{'=' * 80}")
    print("(what fraction of sequences in Rodriguez use at-risk V genes)\n")

    pcp_df = pd.read_csv(PCP_PATH)
    # Filter to depth=2 (MRCA sequences, as used in the entrenchment analysis)
    pcp_depth2 = pcp_df[pcp_df["depth"] == 2].copy()
    print(f"Total depth=2 PCPs in Rodriguez: {len(pcp_depth2)}")

    for v_family in WITHIN_FAMILIES:
        family_pcps = pcp_depth2[pcp_depth2["v_family"] == v_family]
        n_family_pcps = len(family_pcps)
        if n_family_pcps == 0:
            continue

        ent_sites = sorted(entrenched_sites_by_family.get(v_family, set()))
        print(f"\n{v_family} ({n_family_pcps} depth=2 PCPs):")

        for site in ent_sites:
            key = (site, v_family)
            risk_genes = at_risk_genes_per_site.get(key, set())
            if not risk_genes:
                print(f"  Site {site}: 0 at-risk PCPs (0%)")
                continue

            # Match on gene name without allele (*XX) as well as exact match,
            # since partis might report at a different allele resolution.
            risk_genes_base = {g.split("*")[0] for g in risk_genes}

            # Count PCPs using at-risk genes (exact allele match)
            exact_match = family_pcps["v_gene"].isin(risk_genes)
            # Also count by base gene name (in case allele resolution differs)
            base_gene = family_pcps["v_gene"].str.split("*").str[0]
            base_match = base_gene.isin(risk_genes_base)

            n_exact = exact_match.sum()
            n_base = base_match.sum()
            pct_exact = 100 * n_exact / n_family_pcps
            pct_base = 100 * n_base / n_family_pcps

            print(f"  Site {site}: {n_exact} PCPs exact allele match ({pct_exact:.1f}%), "
                  f"{n_base} PCPs base gene match ({pct_base:.1f}%)")
            if risk_genes:
                # Show per-gene breakdown
                for gene in sorted(risk_genes):
                    n_g = (family_pcps["v_gene"] == gene).sum()
                    gene_base = gene.split("*")[0]
                    n_gb = (base_gene == gene_base).sum()
                    if n_g > 0 or n_gb > 0:
                        print(f"    {gene}: {n_g} exact, {n_gb} base-gene PCPs")

    # Save
    output_path = "_output/germline_misassignment_risk.csv"
    rows = []
    for pair in confusable_pairs:
        site_diffs = pair["differing_sites"]
        differing_sites_str = ", ".join(sorted(site_diffs.keys()))
        entrenched_diffs = []
        for fam in [pair["family1"], pair["family2"]]:
            if fam in entrenched_sites_by_family:
                entrenched_diffs.extend(
                    s for s in site_diffs if s in entrenched_sites_by_family[fam]
                )
        rows.append({
            "gene1": pair["gene1"],
            "gene2": pair["gene2"],
            "family1": pair["family1"],
            "family2": pair["family2"],
            "total_nt_distance": pair["total_nt_distance"],
            "differing_sites": differing_sites_str,
            "entrenched_sites_differing": ", ".join(sorted(set(entrenched_diffs))),
        })
    pd.DataFrame(rows).to_csv(output_path, index=False)
    print(f"\nFull results saved to {output_path}")


if __name__ == "__main__":
    main()
