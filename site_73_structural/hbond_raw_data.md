# Hydrogen-bond data: Chothia sites 73–75 in IGHV1 and IGHV3 Fab structures

Primary data record for the FR3 hydrogen-bond analysis: every H-bond returned by
ChimeraX for sites 73–75 in each analyzed structure, with donor and acceptor
atoms and donor–acceptor distances, plus density sampling for the two structures
displayed with electron density maps.

This file records measurements only. Summary statistics, cross-family
comparison, and interpretation are in the manuscript (Methods, Supplementary
Table of sites 73–75, and the figure captions), so that each derived number has
a single authoritative source.

---

## Provenance

- **Software:** UCSF ChimeraX 1.11.1 (2026-01-23).
- **Command** (per structure, with chain and residue range substituted):

  ```
  hbonds #1/<chain>:<sites 73-75> restrict any intraRes false reveal false log true
  ```

  `restrict any` reports every H-bond in which sites 73–75 participate as donor
  or acceptor, against any partner in the structure — framework, CDRs, antigen,
  or water. `intraRes false` suppresses within-residue bonds.
- **Geometry:** ChimeraX default criteria, i.e. idealized per-atom-type geometry
  relaxed by 0.4 Å and 20°, giving an effective absolute D···A cutoff of ~3.5 Å
  for N/O pairs. Distances listed are D···A. Bonds at or beyond ~3.5 Å are
  flagged inline as borderline: they are reported because the cutoff is
  per-atom-type rather than a flat threshold, and would drop out under stricter
  parameters.
- **Chain selection:** one heavy chain per structure. All heavy-chain copies
  within each entry are sequence-identical, so the choice is for
  non-redundancy, not selection of a preferred conformer. The one exception is
  documented under 6VY4.
- **Classification:** backbone (bb) = amide N and carbonyl O. All other donor
  and acceptor atoms are treated as sidechain (sc).
- **Reciprocal pairs:** Thr–Thr and Thr–water hydroxyl contacts are returned
  twice by ChimeraX, once in each direction, because either partner can donate.
  Both directions are listed as returned, with the duplicate marked, so counts
  here reconcile with the ChimeraX log. They are one physical contact.

### Residue numbering

Site labels are **Chothia**; the Donor and Acceptor columns give the **PDB**
atom specification verbatim, so every row is traceable to the raw output.

PDB→Chothia offsets at FR3 (sites 71–77):

| PDB  | Chain | FR3 offset | Chothia 73 = PDB |
|------|-------|-----------|------------------|
| 4H8W | H     | 0         | 73               |
| 3BN9 | D     | 0         | 73               |
| 6ULE | A     | 0         | 73               |
| 8G3Z | E     | +1        | 74               |
| 7X29 | F     | +1        | 74               |
| 6VY4 | H     | +1        | 74               |
| 6PPG | B     | +1        | 74               |
| 2NY1 | D     | +3001     | 3074             |

**The offset is not uniform along the chain.** In the +1 and +3001 structures
the shift arises from a CDR insertion numbered sequentially rather than with an
insertion code, so positions N-terminal to that insertion carry a smaller
offset. Two consequences matter here:

- **K23**, the lysine contacting T75 in some IGHV1 structures, is at PDB 23 in
  8G3Z (offset 0 at that position) and PDB 3023 in 2NY1 (offset +3000).
  Applying the FR3 offset to position 23 selects the wrong residue.
- **CDR-H2 position 52A** is numbered `52A` in 4H8W, 3BN9 and 6ULE, which use
  insertion codes, but `53` in 6VY4 and 6PPG, which do not: 6PPG chain B runs
  52 ASN / 53 TRP / 54 SER and 6VY4 chain H runs 52 THR / 53 PRO / 54 ILE. The
  CDR-H2 acceptor is at Chothia **52A in all four** IGHV3 structures; only the
  residue identity varies (His / Asn / Gly / Trp).

---

## IGHV1 — sites K73, S74, T75

### 2NY1 — 17b Fab + HIV-1 gp120 · chain D · X-ray 1.99 Å

Query: `#1/D:3074-3076`. Motif ADKST at PDB 3072–3076.

**Protein–protein**

| Donor (PDB) | Acceptor (PDB) | Chothia | Type | D···A (Å) |
|---|---|---|---|---|
| /D LYS 3023 NZ  | /D THR 3076 O   | K23 sc → T75 bb  | sc→bb | 3.011 |
| /D SER 3075 N   | /D ASP 3073 OD2 | S74 bb → D72 sc  | bb→sc | 3.296 |
| /D SER 3075 OG  | /D ASP 3073 OD2 | S74 sc → D72 sc  | sc→sc | 2.759 |
| /D THR 3076 N   | /D ASP 3073 O   | T75 bb → D72 bb  | bb→bb | 3.320 |
| /D THR 3076 OG1 | /D THR 3078 OG1 | T75 sc → T77 sc  | sc→sc | 3.084 |
| /D THR 3078 OG1 | /D THR 3076 O   | T77 sc → T75 bb  | sc→bb | 2.828 |
| /D THR 3078 OG1 | /D THR 3076 OG1 | T77 sc → T75 sc  | sc→sc | 3.084 *(reciprocal of row 5)* |

**Water-mediated**

| Donor (PDB) | Acceptor (PDB) | Chothia | D···A (Å) |
|---|---|---|---|
| /D HOH 200 O    | /D LYS 3074 O   | water → K73 bb | 2.328 |
| /D LYS 3074 N   | /D HOH 378 O    | K73 bb → water | 2.920 |
| /D HOH 116 O    | /D THR 3076 OG1 | water → T75 sc | 2.750 |
| /D THR 3076 OG1 | /D HOH 116 O    | T75 sc → water | 2.750 *(reciprocal)* |

**Notes.** K73 makes no protein H-bond of any kind. Its only contacts are two
waters, both to backbone atoms (N, O); the Lys NZ makes no H-bond at all.
T75 sc contacts T77 but not D72 here — the T75–D72 contact is
backbone-to-backbone.

#### Density at K73

2mFo–DFc map (`open 2ny1 from eds`), map σ = 0.2369 in map units. Value
interpolated at each atom centre, expressed in σ:

| Atom | σ |   | Atom | σ |
|---|---|---|---|---|
| N  | 3.12 | | CG | 1.76 |
| CA | 3.27 | | CD | 1.14 |
| C  | 3.83 | | CE | **0.97** |
| O  | 1.72 | | NZ | **0.86** |
| CB | 1.90 | |    |      |

All nine atoms of the residue are modeled. At a 1.0σ contour, Cε and Nζ fall
below threshold. Framework backbone atoms in this window are ≥2.2σ and T77 is
well ordered throughout (OG1 = 2.81σ). Every displayed atom lies inside the
deposited map box, and the Cδ → Cε → Nζ falloff (1.14 → 0.97 → 0.86) is a
genuine weakening of density rather than a map-extent effect (see the note
under 4H8W).

---

### 6VY4 — HENV-32 Fab + Hendra receptor binding protein · chain H · X-ray 2.00 Å

Query: `#1/H:74-76`.

**Chain choice.** Chain H is used rather than chain C. Chain C has a truncated
K73 sidechain — modeled only to Cβ, with no CG/CD/CE/NZ — which makes any K73
sidechain H-bond undetectable there. Chain H is sequence-identical and has the
complete Lys sidechain. This is a local density artifact in chain C, not a
structural difference.

**Protein–protein**

| Donor (PDB) | Acceptor (PDB) | Chothia | Type | D···A (Å) |
|---|---|---|---|---|
| /H LYS 74 NZ   | /H PRO 53 O   | K73 sc → CDR-H2 Pro52A bb | sc→bb | 3.474 |
| /H THR 76 N    | /H ASP 73 O   | T75 bb → D72 bb | bb→bb | 3.168 |
| /H THR 76 OG1  | /H ASP 73 O   | T75 sc → D72 bb | sc→bb | 3.456 |
| /H THR 78 OG1  | /H THR 76 O   | T77 sc → T75 bb | sc→bb | 3.134 |
| /H THR 78 OG1  | /H THR 76 OG1 | T77 sc → T75 sc | sc→sc | 3.136 |

No water-mediated bonds returned.

**Notes.** S74 makes no H-bonds in chain H; the S74 rotamer here does not reach
D72. Chain C has S74 OG → D72 OD1 at 2.39 Å, a rotamer difference between
copies that is independent of the K73 truncation. No K23 contact in this
structure. The K73 (3.474 Å) and T75 sc → D72 bb (3.456 Å) contacts are both
close to the effective cutoff.

---

### 8G3Z — FNI17 Fab + influenza B neuraminidase · chain E · cryo-EM 2.30 Å

Query: `#1/E:74-76`.

**Protein–protein**

| Donor (PDB) | Acceptor (PDB) | Chothia | Type | D···A (Å) |
|---|---|---|---|---|
| /E LYS 23 NZ  | /E THR 76 O   | K23 sc → T75 bb | sc→bb | 3.007 |
| /E SER 75 N   | /E ASP 73 OD1 | S74 bb → D72 sc | bb→sc | 2.930 |
| /E SER 75 OG  | /E ASP 73 OD1 | S74 sc → D72 sc | sc→sc | 2.551 |
| /E SER 75 OG  | /E ASP 73 OD2 | S74 sc → D72 sc | sc→sc | 3.563 *(borderline, >3.5 Å)* |
| /E THR 76 N   | /E ASP 73 OD1 | T75 bb → D72 sc | bb→sc | 3.392 |
| /E THR 76 OG1 | /E ASP 73 OD2 | T75 sc → D72 sc | sc→sc | 2.626 |
| /E THR 76 OG1 | /E THR 78 OG1 | T75 sc → T77 sc | sc→sc | 2.762 |
| /E THR 78 OG1 | /E THR 76 OG1 | T77 sc → T75 sc | sc→sc | 2.762 *(reciprocal of row 7)* |

No water-mediated bonds; no ordered solvent is modeled in this cryo-EM entry.

**Notes.** K73 makes no H-bonds. There is no T77 → T75 backbone contact here,
unlike the other three IGHV1 structures. K23 is at PDB 23 (offset 0 at that
position).

---

### 7X29 — S41 Fab + MERS-CoV spike · chain F · cryo-EM 2.49 Å

Query: `#1/F:74-76`.

**Protein–protein**

| Donor (PDB) | Acceptor (PDB) | Chothia | Type | D···A (Å) |
|---|---|---|---|---|
| /F LYS 74 NZ  | /C ASP 510 O  | K73 sc → antigen bb | sc→bb | 3.475 |
| /F THR 76 N   | /F ASP 73 OD2 | T75 bb → D72 sc | bb→sc | 3.279 |
| /F THR 76 OG1 | /F ASP 73 OD2 | T75 sc → D72 sc | sc→sc | 2.911 |
| /F SER 77 OG  | /F THR 76 O   | S76 sc → T75 bb | sc→bb | 2.699 |
| /F THR 78 OG1 | /F THR 76 O   | T77 sc → T75 bb | sc→bb | 3.407 |

No water-mediated bonds returned.

**Notes.** The K73 partner is Asp510 of the MERS-CoV spike glycoprotein
(chain C) — an intermolecular contact with antigen, not framework. At 3.475 Å
it is close to the effective cutoff. S74 makes no H-bonds. T75 sc → T77 sc is
absent here, the only IGHV1 structure where that is the case.

---

## IGHV3 — sites N73, S74/A74, K75

### 4H8W — N5-i5 Fab + HIV-1 gp120 · chain H · X-ray 1.85 Å

Query: `#1/H:73-75`. Motif RDNSK.

**Protein–protein**

| Donor (PDB) | Acceptor (PDB) | Chothia | Type | D···A (Å) |
|---|---|---|---|---|
| /H ARG 71 NE  | /H ASN 73 OD1 | R71 sc → N73 sc | sc→sc | 2.843 |
| /H ARG 71 NH2 | /H ASN 73 OD1 | R71 sc → N73 sc | sc→sc | 3.496 *(borderline)* |
| /H ASN 73 ND2 | /H ASN 52A O  | N73 sc → CDR-H2 Asn52A bb | sc→bb | 3.036 |
| /H SER 74 N   | /H ASP 72 OD1 | S74 bb → D72 sc | bb→sc | 2.924 |
| /H SER 74 OG  | /H ASP 72 OD1 | S74 sc → D72 sc | sc→sc | 2.547 *(alt loc B, occ. 0.58)* |
| /H LYS 75 N   | /H ASP 72 O   | K75 bb → D72 bb | bb→bb | 3.173 |
| /H LYS 75 N   | /H ASP 72 OD1 | K75 bb → D72 sc | bb→sc | 3.254 |
| /H THR 77 OG1 | /H LYS 75 O   | T77 sc → K75 bb | sc→bb | 3.130 |

**Water-mediated** (5)

| Donor (PDB) | Acceptor (PDB) | Chothia | D···A (Å) |
|---|---|---|---|
| /H ASN 73 N   | /H HOH 350 O | N73 bb → water | 2.919 |
| /H ASN 73 ND2 | /H HOH 447 O | N73 sc → water | 3.164 |
| /H ASN 73 ND2 | /H HOH 475 O | N73 sc → water | 2.764 |
| /H HOH 436 O  | /H ASN 73 O  | water → N73 bb | 3.174 |
| /H HOH 461 O  | /H ASN 73 O  | water → N73 bb | 3.380 |

**Notes.** K75 NZ makes no H-bond; all K75 contacts are backbone. The S74 OG
bond derives from alternate location B (occupancy 0.58).

#### Density at R71 and N73, and a note on map coverage

2mFo–DFc map (`open 4h8w from eds`), deposited map σ = 0.2210.

Sampling the deposited map directly at atom centres:

| Residue | Atom | σ | | Residue | Atom | σ |
|---|---|---|---|---|---|---|
| R71 | N   | 4.14 | | N73 | N   | 4.06 |
| R71 | CA  | 3.46 | | N73 | CA  | 3.76 |
| R71 | C   | 4.59 | | N73 | C   | 3.35 |
| R71 | O   | 3.55 | | N73 | O   | 2.19 |
| R71 | CB  | 3.01 | | N73 | CB  | 2.90 |
| R71 | CG  | 3.52 | | N73 | CG  | 3.88 |
| R71 | CD  | 0.00 | | N73 | OD1 | 3.91 |
| R71 | NE  | 0.00 | | N73 | ND2 | 2.98 |
| R71 | CZ  | 0.00 | | | | |
| R71 | NH1 | 0.00 | | | | |
| R71 | NH2 | 0.00 | | | | |

The five zero values are a map-extent effect, not weak density. The deposited
map covers exactly one unit cell — grid 100 × 120 × 150 at steps
0.9225 / 0.8659 / 0.8977 Å, giving a = 92.251, b = 103.909, c = 134.654 Å with
all cell angles 90°, spanning x ∈ [0, 91.33]. R71 Cδ, Nε, Cζ, Nη1 and Nη2 sit
at x = −1.3 to −4.0 Å, so the guanidinium protrudes past the x = 0 face into
the adjacent unit cell, where the deposited grid holds no samples. Interpolation
outside the grid returns exactly zero.

The density for those atoms is present in the deposited file, one lattice
translation away. Sampling at x + a reproduces it:

| Atom | at deposited coords | at x + a (92.251 Å) |
|---|---|---|
| R71 CD  | 0.00 | 2.54 |
| R71 NE  | 0.00 | 4.74 |
| R71 CZ  | 0.00 | 5.04 |
| R71 NH1 | 0.00 | 2.90 |
| R71 NH2 | 0.00 | 3.50 |

The converse also holds — atoms inside the box return 0.00 when shifted by +a —
confirming the map spans one cell with no padding. This is plain periodic
tiling, not a rotated symmetry mate.

ChimeraX tiles the map over a chosen region natively, using the cell and space
group in the map header:

```
volume cover #2 atomBox #1/H:71-77 pad 5
```

This reports 4 symmetry operations, 0 points not covered by any symmetry, and a
maximum discrepancy of 3.3 × 10⁻⁶ where copies overlap. Re-sampling the covered
map reproduces the x + a values above exactly.

**Interpretation for display.** The R71 guanidinium is among the best-ordered
atoms in this loop (Nε 4.74σ, Cζ 5.04σ — stronger than any atom of N73), and it
makes the 2.84 Å bond to N73 OD1. Any apparent absence of density around it in
a figure rendered from the deposited map reflects the cell boundary, not
disorder.

**If rendering with `volume cover`:** the covered map has its own σ (0.2765),
larger than the deposited map's (0.2210), because the covered box is small and
protein-centred. Contouring it with `sdLevel 1.0` would therefore sit at
~1.25σ in deposited-map terms. To contour at 1.0σ of the deposited map, set the
level absolutely:

```
volume #3 style surface level 0.221
```

---

### 3BN9 — E2 Fab + MT-SP1 · chain D · X-ray 2.17 Å

Query: `#1/D:73-75`. Motif RDNSK.

**Protein–protein**

| Donor (PDB) | Acceptor (PDB) | Chothia | Type | D···A (Å) |
|---|---|---|---|---|
| /D ARG 71 NE  | /D ASN 73 OD1 | R71 sc → N73 sc | sc→sc | 2.910 |
| /D ARG 71 NH1 | /D ASN 73 OD1 | R71 sc → N73 sc | sc→sc | 3.513 *(borderline, >3.5 Å)* |
| /D ASN 73 ND2 | /D GLY 52A O  | N73 sc → CDR-H2 Gly52A bb | sc→bb | 2.954 |
| /D SER 74 N   | /D ASP 72 OD1 | S74 bb → D72 sc | bb→sc | 3.115 |
| /D SER 74 OG  | /D ASP 72 OD1 | S74 sc → D72 sc | sc→sc | 2.606 |
| /D LYS 75 N   | /D ASP 72 O   | K75 bb → D72 bb | bb→bb | 3.078 |
| /D LYS 75 N   | /D ASP 72 OD1 | K75 bb → D72 sc | bb→sc | 3.355 |
| /D THR 77 OG1 | /D LYS 75 O   | T77 sc → K75 bb | sc→bb | 3.021 |

**Water-mediated** (1)

| Donor (PDB) | Acceptor (PDB) | Chothia | D···A (Å) |
|---|---|---|---|
| /D ASN 73 N | /D HOH 272 O | N73 bb → water | 3.011 |

**Notes.** K75 NZ makes no H-bond.

---

### 6ULE — 2541 Fab + *P. falciparum* circumsporozoite protein · chain A · X-ray 2.55 Å

Query: `#1/A:73-75`. Motif RDNSK.

**Protein–protein**

| Donor (PDB) | Acceptor (PDB) | Chothia | Type | D···A (Å) |
|---|---|---|---|---|
| /A ARG 71 NH2 | /A ASN 73 OD1 | R71 sc → N73 sc | sc→sc | 3.133 |
| /A ASN 73 ND2 | /A HIS 52A O  | N73 sc → CDR-H2 His52A bb | sc→bb | 3.209 |
| /A LYS 75 N   | /A ASP 72 O   | K75 bb → D72 bb | bb→bb | 3.274 |

**Water-mediated** (1)

| Donor (PDB) | Acceptor (PDB) | Chothia | D···A (Å) |
|---|---|---|---|
| /A ASN 73 ND2 | /A HOH 306 O | N73 sc → water | 2.877 |

**Notes.** The sparsest of the four IGHV3 structures. Only one of the two R71
guanidinium contacts is present (NH2 only; no NE). S74 makes no H-bonds. K75
has only the single backbone contact to D72 — no K75 bb → D72 sc, and no
T77/S77 → K75 contact. K75 NZ makes no H-bond.

---

### 6PPG — MCAF5352A Fab + IL-17F · chain B · X-ray 2.75 Å

Query: `#1/B:74-76`. Motif **RDNAK** — site 74 is Ala, which is the germline
residue for this V gene, so the unmutated-germline selection criterion is
satisfied.

**Protein–protein**

| Donor (PDB) | Acceptor (PDB) | Chothia | Type | D···A (Å) |
|---|---|---|---|---|
| /B ARG 72 NE  | /B ASN 74 OD1 | R71 sc → N73 sc | sc→sc | 2.836 |
| /B ARG 72 NH1 | /B ASN 74 OD1 | R71 sc → N73 sc | sc→sc | 3.407 |
| /B ASN 74 ND2 | /B TRP 53 O   | N73 sc → CDR-H2 Trp52A bb | sc→bb | 3.049 |
| /B LYS 76 N   | /B ASP 73 O   | K75 bb → D72 bb | bb→bb | 3.191 |
| /B LYS 76 N   | /B ASP 73 OD1 | K75 bb → D72 sc | bb→sc | 2.706 |
| /B SER 78 OG  | /B LYS 76 O   | S77 sc → K75 bb | sc→bb | 2.775 |

**Water-mediated** (1)

| Donor (PDB) | Acceptor (PDB) | Chothia | D···A (Å) |
|---|---|---|---|
| /B ASN 74 N | /B HOH 405 O | N73 bb → water | 2.724 |

**Notes.** A74 has no sidechain hydroxyl and makes no H-bonds. Position 77 is
Ser here rather than Thr, and still donates to the K75 backbone carbonyl. K75
NZ makes no H-bond. The CDR-H2 acceptor is Trp at PDB 53 = Chothia 52A; this
entry does not use insertion codes.
