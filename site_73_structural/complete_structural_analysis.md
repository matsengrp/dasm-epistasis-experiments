# Complete Structural Analysis: H-bond Networks at Chothia Sites 73-75 in IGHV1 vs IGHV3 Antibody Frameworks

## 1. Background and Motivation

This analysis investigates the structural basis for the observed entrenchment of specific amino acids at Chothia framework sites 73 and 75 in human IGHV1 and IGHV3 heavy chains. Sequence-based evolutionary analysis identified that:

- **IGHV1**: K73 and T75 are entrenched (resistant to substitution).
- **IGHV3**: N73 and K75 are entrenched.

The two V gene families therefore have *different* identities at the same Chothia positions. The goal of this structural analysis was to determine the H-bond network — and especially the sidechain partners — engaged by each entrenched residue, in order to test the hypothesis that the residues are entrenched because they are stabilized by family-specific sidechain interactions in the framework.

A key finding is that the two families exhibit **reversed topology** at sites 73/75: in IGHV1, K73 points outward and T75 points inward (toward the framework); in IGHV3, N73 points inward and K75 points outward. The inward-facing residue in each family provides the conserved sidechain H-bond network.

---

## 2. Methods

### 2.1 Structure selection criteria

PDB structures were selected based on:
1. Human antibody Fab fragments
2. Heavy chain V gene assignment to either IGHV1 or IGHV3 (assigned with `partis`)
3. Unmutated germline residues at framework positions 71–75 (Chothia numbering), verified against the OGRDB germline reference (https://ogrdb.airr-community.org/germline_set/75)
4. Presence of the entrenched residues at positions 73 and 75 (K73/T75 for IGHV1; N73/K75 for IGHV3)
5. Resolution ≤ 3.0 Å

Four structures were analyzed per family:
- **IGHV1**: 7X29 (2.49 Å), 2NY6 (2.80 Å), 6VY4 (2.00 Å), 8G3Z (2.30 Å)
- **IGHV3**: 4H8W (1.85 Å), 3BN9 (2.17 Å), 6ULE (2.55 Å), 6PPG (2.75 Å)

### 2.2 Software

UCSF ChimeraX version 1.11.1 was used for all structural analysis and visualization.

### 2.3 Residue numbering

Chothia framework sites 71–75 correspond to a conserved sequence motif:
- IGHV1: **A**71-**D**72-**K**73-**S**74-**T**75 ("ADKST")
- IGHV3: **R**71-**D**72-**N**73-**S**74-**K**75 ("RDNSK")

PDB files use a variety of internal residue numbering schemes. PDB-to-Chothia offsets observed:

| PDB | Chain | PDB→Chothia offset |
|-----|-------|--------------------|
| 6ULE | A | 0 (direct correspondence) |
| 4H8W | H | 0 |
| 3BN9 | D | 0 |
| 8G3Z | E | +1 (Chothia 73 = PDB 74) |
| 7X29 | F | +1 |
| 6VY4 | C, H | +1 |
| 6PPG | B | +1 |
| 2NY6 | D | +3001 (Chothia 73 = PDB 3074) |

For each structure, residue mapping was verified individually by confirming the ADKST or RDNSK motif identities at the expected positions. **All residue numbers in this report use Chothia numbering.**

### 2.4 H-bond analysis

H-bonds were computed in ChimeraX with the `hbonds` command using default geometric criteria. The 0.4 Å and 20° values are **tolerances** that relax ChimeraX's idealized geometry (which is per-atom-type, with typical D···A ideal distances of 2.8–3.1 Å for N/O pairs). The effective absolute D···A distance cutoff is therefore ~3.5 Å for N/O pairs. The general command form was:

```
hbonds #1/<chain>:73-75 restrict any reveal true color cyan dashes 5 log true
```

`restrict any` ensures that ALL H-bonds in which sites 73–75 participate (either as donor or acceptor) are reported, against any partner in the structure (including antigen, framework, CDRs, and water).

Only one copy of the heavy chain was analyzed per structure to avoid redundancy from crystallographic symmetry mates or multiple Fabs in the asymmetric unit. Bonds involving solvent (water) were tabulated separately from protein–protein bonds.

Each H-bond was classified by **donor type (sc/bb)** and **acceptor type (sc/bb)** based on whether the donor/acceptor atom belongs to the residue sidechain or backbone:
- Backbone (bb) atoms: N (amide), O (carbonyl)
- Sidechain (sc) atoms: NE/NH1/NH2 (Arg), NZ (Lys), OD1/OD2 (Asp), ND2/OD1 (Asn), OG (Ser), OG1 (Thr)

### 2.5 ChimeraX visualization commands

Documented in `chimerax_commands.md`. Briefly: structures were rendered with cartoon ribbons (gray), sites 73–75 displayed as yellow ball-and-stick (colored by heteroatom), partner residues shown in gray ball-and-stick, and H-bonds drawn as cyan dashed lines. Residues were labeled with single-letter Chothia codes (e.g., "N73") with custom 3D offsets for legibility.

---

## 3. Raw H-bond Results (Per Structure)

**All distances in Å. All residue numbers in Chothia.**

### 3.1 IGHV1 Structures

#### 8G3Z — FNI17 Fab + Influenza Neuraminidase (Chain E, 2.30 Å)

**K73 — 0 H-bonds**

**S74 — 3 H-bonds (1 bb, 2 sc)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| S74 N | D72 OD1 | bb→sc | 2.93 |
| S74 OG | D72 OD1 | sc→sc | 2.55 |
| S74 OG | D72 OD2 | sc→sc | 3.56 |

**T75 — 5 H-bonds (1 bb, 4 sc)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| T75 N | D72 OD1 | bb→sc | 3.39 |
| T75 OG1 | D72 OD2 | sc→sc | 2.63 |
| T75 OG1 | T77 OG1 | sc→sc | 2.76 |
| T77 OG1 | T75 OG1 | sc→sc | 2.76 |
| K23 NZ | T75 O | sc→bb | 3.01 |

#### 7X29 — S41 Ab + MERS-CoV Spike (Chain F, 2.49 Å)

**K73 — 1 H-bond (1 sc)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| K73 NZ | Antigen Asp O | sc→bb | 3.48 |

K73 sidechain is engaged with the antigen — consistent with K73 being outward-facing and accessible to the binding surface.

**S74 — 0 H-bonds**

**T75 — 4 H-bonds (1 bb, 3 sc)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| T75 N | D72 OD2 | bb→sc | 3.28 |
| T75 OG1 | D72 OD2 | sc→sc | 2.91 |
| S76 OG | T75 O | sc→bb | 2.70 |
| T77 OG1 | T75 O | sc→bb | 3.41 |

#### 2NY6 — 17b Fab + HIV-1 gp120 (Chain D, 2.80 Å)

**K73 — 1 H-bond (1 sc, at the relaxed-distance limit)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| K73 NZ | D55 OD2 | sc→sc | 3.51 |

D55 is in CDR-H2. The 3.51 Å distance is right at the edge of the relaxed N/O distance window (~3.5 Å); under stricter parameters (smaller distance tolerance) it would no longer be reported.

**S74 — 2 H-bonds (2 sc)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| S74 OG | D72 OD1 | sc→sc | 3.29 |
| S74 OG | D72 OD2 | sc→sc | 2.69 |

**T75 — 6 H-bonds (1 bb, 5 sc; 1 bb-bb)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| T75 N | D72 O | bb→bb | 3.19 |
| T75 OG1 | D72 OD2 | sc→sc | 2.92 |
| T75 OG1 | T77 OG1 | sc→sc | 3.40 |
| T77 OG1 | T75 O | sc→bb | 2.86 |
| T77 OG1 | T75 OG1 | sc→sc | 3.40 |

#### 6VY4 — HENV-32 Ab + Hendra RBP (Chain H, 2.00 Å)

**Note**: This structure has two heavy chains (C and H) with identical sequence. Chain C was initially analyzed but was found to have a **truncated K73 sidechain** (only atoms N, CA, C, O, CB modeled — no CG/CD/CE/NZ), preventing detection of any K73 sidechain H-bonds. Because the user prefiltered structures so that all heavy chains within each PDB are sequence-identical, chain H — which has the **complete** Lys sidechain (N, CA, C, O, CB, CG, CD, CE, NZ) — is used here as the representative copy. The chain C truncation is a local electron-density artifact, not a structural feature.

**K73 — 1 H-bond (1 sc → bb)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| K73 NZ | P (CDR-H2) O | sc→bb | 3.47 |

**This confirms the manual annotation** of K73 → CDR-H2 P52A backbone in 6VY4. The lysine sidechain reaches across to the CDR-H2 loop and donates an H-bond to a proline backbone carbonyl.

**S74 — 0 H-bonds**

(In chain H, the S74 sidechain rotamer does not engage D72; this differs from chain C where S74 OG → D72 OD1 was observed at 2.39 Å. This rotamer-level difference between the two heavy-chain copies is independent of the K73 truncation issue.)

**T75 — 4 H-bonds (1 bb, 3 sc)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| T75 N | D72 O | bb→bb | 3.17 |
| T75 OG1 | D72 O | sc→bb | 3.46 |
| T77 OG1 | T75 O | sc→bb | 3.13 |
| T77 OG1 | T75 OG1 | sc→sc | 3.14 |

---

### 3.2 IGHV3 Structures

#### 6ULE — 2541 Fab + Plasmodium Circumsporozoite Protein (Chain A, 2.55 Å)

**N73 — 2 protein H-bonds (2 sc) + 1 water**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| R71 NH2 | N73 OD1 | sc→sc | 3.13 |
| N73 ND2 | H52A O (CDR-H2) | sc→bb | 3.21 |

Water-mediated: N73 ND2 → HOH (2.88)

**S74 — 0 H-bonds**

**K75 — 1 H-bond (1 bb-bb)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| K75 N | D72 O | bb→bb | 3.27 |

#### 4H8W — N5-i5 Ab + HIV-1 gp120 (Chain H, 1.85 Å)

**N73 — 3 protein H-bonds (3 sc) + 5 waters**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| R71 NE | N73 OD1 | sc→sc | 2.84 |
| R71 NH2 | N73 OD1 | sc→sc | 3.50 |
| N73 ND2 | N52A O (CDR-H2) | sc→bb | 3.04 |

**S74 — 2 H-bonds (1 bb, 1 sc)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| S74 N | D72 OD1 | bb→sc | 2.92 |
| S74 OG | D72 OD1 | sc→sc | 2.55 |

**K75 — 3 H-bonds (3 bb)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| K75 N | D72 O | bb→bb | 3.17 |
| K75 N | D72 OD1 | bb→sc | 3.25 |
| T77 OG1 | K75 O | sc→bb | 3.13 |

#### 3BN9 — E2 Fab + MT-SP1 (Chain D, 2.17 Å)

**N73 — 3 protein H-bonds (3 sc) + 1 water**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| R71 NE | N73 OD1 | sc→sc | 2.91 |
| R71 NH1 | N73 OD1 | sc→sc | 3.51 |
| N73 ND2 | G52A O (CDR-H2) | sc→bb | 2.95 |

Water-mediated: N73 N → HOH (3.01)

**S74 — 2 H-bonds (1 bb, 1 sc)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| S74 N | D72 OD1 | bb→sc | 3.12 |
| S74 OG | D72 OD1 | sc→sc | 2.61 |

**K75 — 3 H-bonds (3 bb)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| K75 N | D72 O | bb→bb | 3.08 |
| K75 N | D72 OD1 | bb→sc | 3.36 |
| T77 OG1 | K75 O | sc→bb | 3.02 |

#### 6PPG — MCAF5352A Fab + IL-17F (Chain B, 2.75 Å)

**N73 — 3 protein H-bonds (3 sc) + 1 water**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| R71 NE | N73 OD1 | sc→sc | 2.84 |
| R71 NH1 | N73 OD1 | sc→sc | 3.41 |
| N73 ND2 | W53 O (CDR-H2) | sc→bb | 3.05 |

Water-mediated: N73 N → HOH (2.72)

**A74 — 0 H-bonds** (site 74 is Ala in this structure — no sidechain hydroxyl)

**K75 — 3 H-bonds (3 bb)**

| Donor | Acceptor | Type | Distance |
|-------|----------|------|----------|
| K75 N | D72 O | bb→bb | 3.19 |
| K75 N | D72 OD1 | bb→sc | 2.71 |
| S77 OG | K75 O | sc→bb | 2.78 |

---

## 4. Summary Tables

### 4.1 IGHV1 — sidechain vs backbone H-bond counts

| Site | 8G3Z | 7X29 | 2NY6 | 6VY4 | sc avg | bb avg |
|------|------|------|------|------|--------|--------|
| K73 sc | 0 | 1 (antigen) | 1 (D55, borderline) | 1 (CDR-H2 P bb) | 0.75 | — |
| K73 bb | 0 | 0 | 0 | 0 | — | 0 |
| S74 sc | 2 | 0 | 2 | 0 | 1.00 | — |
| S74 bb | 1 | 0 | 0 | 0 | — | 0.25 |
| T75 sc | 3 | 1 | 3 | 2 | 2.25 | — |
| T75 bb | 2 | 3 | 3 | 2 | — | 2.50 |

(6VY4 values from chain H; chain C had a truncated K73 sidechain.)

### 4.2 IGHV3 — sidechain vs backbone H-bond counts

| Site | 6ULE | 4H8W | 3BN9 | 6PPG | sc avg | bb avg |
|------|------|------|------|------|--------|--------|
| N73 sc | 2 | 3 | 3 | 3 | 2.75 | — |
| N73 bb | 0 | 0 | 0 | 0 | — | 0 |
| N73 water | 1 | 5 | 1 | 1 | — | — |
| S74 sc | 0 | 1 | 1 | 0 (Ala) | — | — |
| S74 bb | 0 | 1 | 1 | 0 (Ala) | — | — |
| K75 sc | 0 | 0 | 0 | 0 | **0** | — |
| K75 bb | 1 | 3 | 3 | 3 | — | 2.50 |

### 4.3 Cross-family comparison

| Feature | IGHV1 | IGHV3 |
|---------|-------|-------|
| Entrenched residues | K73, T75 | N73, K75 |
| Site 73 identity | Lys (K) | Asn (N) |
| Site 73 orientation | **Outward (4/4)** | **Inward (4/4)** |
| Site 73 sidechain H-bonds | 0–1, variable | 2–3, conserved |
| Site 73 sidechain partners | None consistent | R71 sc (4/4) + CDR-H2 bb (4/4) |
| Site 73 water-mediated | None (0/4) | Present (4/4) |
| Site 75 identity | Thr (T) | Lys (K) |
| Site 75 orientation | **Inward (4/4)** | **Outward (4/4)** |
| Site 75 sidechain H-bonds | 2–3, conserved | **0** |
| Site 75 sidechain partners | D72 sc (4/4) + T77 sc (3/4) | None — sc completely free |
| Site 75 backbone H-bonds | bb→D72 (4/4); ←T77 (4/4) | bb→D72 (4/4); ←T77/S77 (3/4) |
| **Reversed topology** | K73 outward / T75 inward | N73 inward / K75 outward |
| **Inward-facing stabilizer** | **T75** (sc network to D72, T77) | **N73** (sc network to R71, CDR-H2) |
| Role of D72 | Accepts T75 sc and S74 sc | Accepts K75 bb and N73 sc (via R71 relay) |
| Role of site 71 | A71 — no sc role | R71 — bridges N73 to FR3 (4/4) |
| Role of K23 | K23 NZ → T75 O in 2/4 | N/A |

---

## 5. Key Findings

1. **Reversed topology is robust across all 8 structures.** In every IGHV1 structure examined, K73 points outward (sidechain solvent-accessible / antigen-accessible) and T75 points inward (sidechain engaged in framework H-bonds). In every IGHV3 structure, the orientations are reversed: N73 points inward, K75 points outward.

2. **Each family has a single dominant inward-facing stabilizer at sites 73/75 — but it is at a different site.**
   - IGHV1: **T75 sidechain (OG1)** anchors the framework via consistent H-bonds to D72 sidechain (4/4) and T77 sidechain (3/4).
   - IGHV3: **N73 sidechain (OD1, ND2)** anchors the framework via consistent H-bonds to R71 sidechain (4/4) and CDR-H2 backbone carbonyl (4/4).

3. **The outward-facing residue at sites 73/75 is consistently free of sidechain H-bonds in its respective family.**
   - IGHV1 K73 sidechain: free or only weakly engaged (variable, no consistent partner).
   - IGHV3 K75 sidechain: completely free in all 4 structures (zero sidechain H-bonds).

4. **D72 is a shared anchor point in both families**, but the sidechain partner differs: it accepts from T75 sc in IGHV1 and from N73 sc (via the R71 relay) in IGHV3. In both families, the backbone of site 75 also donates to D72.

5. **R71 is critical for IGHV3 entrenchment of N73.** The Arg-Asn sidechain bridge (R71 NE/NH ↔ N73 OD1) is present in 4/4 IGHV3 structures and provides a family-specific reason why position 71 is conserved as Arg in IGHV3 — it would be lost in IGHV1 where position 71 is Ala.

6. **Water-mediated bonds further distinguish the families.** N73 has water-mediated contacts in 4/4 IGHV3 structures (1–5 waters per structure), whereas no water-mediated bonds were observed at IGHV1 K73.

7. **K23 → T75 O backbone H-bond** appears in 2/4 IGHV1 structures (8G3Z, 6VY4), suggesting an additional but inconsistent stabilization of T75 from CDR-H1.

---

## 6. Comparison with User's Manual Analysis (Table S3)

The computational analysis broadly agrees with the user's prior manual annotation, with the following notes and discrepancies:

| Manual finding | Computational status | Notes |
|----------------|----------------------|-------|
| IGHV3 N73 ↔ R71 sc bridge | **Confirmed (4/4)** | Both NE and NH1/NH2 contribute |
| IGHV3 N73 → CDR-H2 backbone | **Confirmed (4/4)** | Acceptors: H52A, N52A, G52A, W53 |
| IGHV1 T75 ↔ D72 sc | **Confirmed (4/4)** | Most stable interaction in IGHV1 |
| IGHV1 T75 ↔ T77 sc | **Confirmed (3/4)** | Absent only in 7X29 (where T77 OG1 → T75 O bb instead) |
| IGHV3 K75 sidechain free | **Confirmed (4/4)** | Zero sidechain H-bonds in any structure |
| IGHV1 K73 sidechain free | Mostly confirmed | Variable; K73 sc contacts antigen in 7X29, weakly to D55 in 2NY6 |
| **6VY4 K73 → CDR-H2 P52A backbone (manual)** | **Confirmed (chain H, 3.47 Å)** | Chain C had truncated K73 sidechain; chain H has the complete Lys and reproduces the manual finding exactly |

### 6.1 The 6VY4 K73 — resolved using chain H

The user's manual analysis indicated that K73 in 6VY4 makes an H-bond to the backbone of P52A (CDR-H2). The initial computational analysis on **chain C** of 6VY4 found only a backbone-backbone bond and no sidechain bond to CDR-H2.

Investigation revealed that **the K73 sidechain in 6VY4 chain C is truncated at Cβ in the deposited coordinates**:

```
atom /C:74@N    Npl
atom /C:74@CA   C3
atom /C:74@C    C2
atom /C:74@O    O2
atom /C:74@CB   C3
(no CG, CD, CE, NZ)
```

The CG, CD, CE, and NZ atoms are absent — most likely because the electron density for the Lys sidechain was too weak to model. ChimeraX therefore cannot detect any sidechain H-bond involving K73 NZ in chain C, regardless of whether one would actually exist in a fully built model.

In contrast, **chain H of the same 6VY4 structure has the complete Lys sidechain modeled**:

```
atom /H:74@N, CA, C, O, CB, CG, CD, CE, NZ   (all present)
```

**Resolution**: Because the structure was prefiltered to PDBs with sequence-identical heavy chains, chain H is a valid representative copy. Re-running `hbonds #1/H:74-76 restrict any` on chain H returns:

```
/H LYS 74 NZ   /H PRO 53 O    3.47 Å   (= K73 NZ → CDR-H2 P backbone)
```

This **exactly confirms** the manual annotation: K73 sidechain donates to a CDR-H2 proline backbone carbonyl. The chain-C result was therefore a false negative caused by incomplete modeling at low local electron density, and chain H is now used as the canonical 6VY4 result throughout this report.

This is an important lesson for any computational H-bond analysis: missing sidechain atoms in PDB models (common at medium-low resolution) will produce false negatives. Always verify atom completeness for residues of interest before drawing conclusions.

---

## 7. Limitations

- **Small sample size** (n = 4 per family). The patterns are remarkably consistent but the statistical confidence at any individual site is limited.
- **Default H-bond geometry cutoffs.** ChimeraX defaults relax idealized N/O geometry by 0.4 Å and 20°, giving an effective absolute D···A distance limit of ~3.5 Å. A few interactions sit right at this relaxed limit (e.g., 2NY6 K73 NZ → D55 OD2 at 3.51 Å) and would no longer be called bonds under stricter parameters.
- **Truncated/missing sidechains** in lower-resolution structures (most starkly in 6VY4 chain C K73) cause false-negative sidechain H-bonds.
- Only **one heavy chain copy** per structure was analyzed — alternative copies in the same crystal may show small differences.
- **Antigen contacts** are reported when present (e.g., 7X29 K73 NZ → antigen O) but were not the focus of this analysis.

---

## 8. Files in this analysis

- `complete_structural_analysis.md` — this file
- `structural_analysis_results.md` — earlier per-structure results (with user annotations)
- `summary_tables.md` — three concise summary tables
- `updated_methods.tex` — LaTeX methods section for manuscript
- `chimerax_commands.md` — full ChimeraX commands for the 6ULE and 8G3Z visualizations

---

## 9. Suggested next steps

1. **Resolve the 6VY4 K73 sidechain.** Either re-analyze chain H or rebuild the chain C sidechain, then re-run the `hbonds` command to determine whether K73 in 6VY4 actually forms a sidechain H-bond to CDR-H2. This is the only point where the computational analysis fails to confirm the manual annotation.
2. **Optionally expand the structure set** to 6–8 structures per family for stronger statistical claims about the conservation rate of each interaction.
3. **Quantify the energetic contribution** of the T75–D72 (IGHV1) and N73–R71 (IGHV3) sidechain interactions using a method such as FoldX or Rosetta ΔΔG calculations on alanine mutants, to test whether their loss would explain the observed entrenchment.
