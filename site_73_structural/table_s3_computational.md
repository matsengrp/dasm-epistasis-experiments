# Table S3 (computational): Structural analysis of hydrogen bonding networks at sites 73–75 in IGHV1 and IGHV3 Fab crystal structures

Each row summarizes the local H-bond interactions observed at the entrenched sites 73–75 for one crystal structure, derived from ChimeraX `hbonds` analysis (default geometry, 0.4 Å / 20°). All residue numbers in Chothia. The reversed topology — K73 outward / T75 inward in IGHV1 versus N73 inward / K75 outward in IGHV3 — is conserved across all 8 examined structures. D72 anchor interactions are reported from the perspective of sites 73–75 (the reciprocal D72 column has been omitted to remove redundancy).

**Notation key.** Each H-bond is annotated by whether the donor and acceptor atoms belong to the residue **sidechain (sc)** or **backbone (bb)**. Backbone atoms are N (amide) and O (carbonyl); all other H-bond donors/acceptors are sidechain. Specific atom names appear in parentheses only when needed to disambiguate which sidechain atom is involved (e.g., "R71 sc (NE + NH2)" indicates two distinct sidechain donors of arginine — the ε nitrogen and a terminal guanidinium nitrogen — both of which are sidechain atoms). Atom name reference: NE/NH1/NH2 = arginine sidechain guanidinium nitrogens; NZ = lysine sidechain ε-amine; OD1/OD2 = aspartate/asparagine sidechain oxygens; ND2 = asparagine sidechain amide nitrogen; OG = serine sidechain hydroxyl; OG1 = threonine sidechain hydroxyl. All of these are sidechain atoms.

## IGHV1 Structures

| PDB  | Res. (Å) | Motif | Site 73 (K) | Site 74 (S) | Site 75 (T) |
|------|----------|-------|-------------|-------------|-------------|
| 6VY4 | 2.00 | ADKST | Outward; sc → CDR-H2 Pro bb (3.47 Å, chain H) | Free (chain H rotamer) | bb → D72 bb; sc → D72 bb; sc ↔ T77 sc; bb ← T77 sc |
| 8G3Z | 2.30 | ADKST | Outward; free | sc + bb → D72 sc | sc + bb → D72 sc; sc ↔ T77 sc; bb ← K23 sc |
| 7X29 | 2.49 | ADKST | Outward; sc → antigen bb | Free | sc + bb → D72 sc; bb ← S76 sc, T77 sc |
| 2NY6 | 2.80 | ADKST | Outward; weak sc → D55 sc (3.51 Å, at the relaxed-distance limit) | sc → D72 sc | bb → D72 bb; sc → D72 sc; sc ↔ T77 sc; bb ← T77 sc |

## IGHV3 Structures

| PDB  | Res. (Å) | Motif | Site 73 (N) | Site 74 (S/A) | Site 75 (K) |
|------|----------|-------|-------------|---------------|-------------|
| 4H8W | 1.85 | RDNSK | Inward; sc ← R71 sc (NE + NH2); sc → CDR-H2 N52A bb | S; sc + bb → D72 sc | Outward; sc free; bb → D72 bb + D72 sc; bb ← T77 sc |
| 3BN9 | 2.17 | RDNSK | Inward; sc ← R71 sc (NE + NH1); sc → CDR-H2 G52A bb; water-mediated | S; sc + bb → D72 sc | Outward; sc free; bb → D72 bb + D72 sc; bb ← T77 sc |
| 6ULE | 2.55 | RDNSK | Inward; sc ← R71 sc (NH2); sc → CDR-H2 H52A bb; water-mediated | S; free | Outward; sc free; bb → D72 bb |
| 6PPG | 2.75 | RDNAK | Inward; sc ← R71 sc (NE + NH1); sc → CDR-H2 W53 bb; water-mediated | A (no sidechain hydroxyl) | Outward; sc free; bb → D72 bb + D72 sc; bb ← S77 sc |

---

### Notes on differences from the manual Table S3

- **6VY4 site 73**: The manual table reports K73 → P52A backbone (CDR-H2). Chain C of 6VY4 has a truncated K73 sidechain (only N, CA, C, O, CB modeled — no CG/CD/CE/NZ), so no sidechain bond was detectable from chain C. Re-running the analysis on **chain H** (which has the complete Lys sidechain and is sequence-identical to chain C) recovers exactly the manual annotation: K73 NZ → CDR-H2 Pro backbone O at 3.47 Å. **The manual finding is confirmed.**
- **2NY6 site 73**: K73 NZ → D55 OD2 contact at 3.51 Å. The 0.4 Å figure in the methods is a *tolerance* on idealized geometry — the effective absolute D···A cutoff for an N/O pair is ~3.5 Å, so this bond sits right at the relaxed limit and would not be reported under stricter parameters. K73 in 2NY6 could reasonably be called "free" in a stricter analysis.
- **6ULE site 75**: The computational analysis finds K75 N → D72 O backbone H-bond only (no bb→sc to D72 OD), making it slightly less anchored than the other three IGHV3 structures.
- **All N73 → CDR-H2 backbone contacts** are confirmed (4/4) at residues H52A, N52A, G52A, W53 — different CDR-H2 residues in different structures, but always to the backbone carbonyl one position into the loop.
- **Site 74 in 6PPG** is Ala (motif RDNAK), not Ser — so no S74-OG → D72 bond is possible.
