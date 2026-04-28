# ChimeraX Visualization Commands

## 6ULE — IGHV3 (2541 Fab + Circumsporozoite Protein)

```
# Load structure
open 6ule

# Set white background
set bgColor white

# Hide all atom representations, show cartoon ribbons
hide #1 atoms
cartoon #1

# Show only the heavy chain (chain A)
hide #1/B,C,D,I cartoons

# Color the heavy chain ribbon gray
color #1/A gray cartoons

# Set simple lighting with no depth cueing
lighting simple
lighting depthCue false

# Add silhouette outlines for clarity
graphics silhouettes true width 1

# Increase rendering quality
graphics quality 4

# Show Chothia sites 73-75 (N73, S74, K75) as yellow ball-and-stick
show #1/A:73-75 atoms
style #1/A:73-75 ball
color #1/A:73-75 yellow atoms
color #1/A:73-75 byhet atoms          # color non-carbon atoms by element

# Show Chothia sites 71-72 (R71, D72) as gray ball-and-stick
show #1/A:71-72 atoms
style #1/A:71-72 ball
color #1/A:71-72 gray atoms
color #1/A:71-72 byhet atoms

# Label residues 71-75 with single-letter codes (Chothia numbering)
label #1/A:71 residues text "R71" color black height 0.7 font Arial
label #1/A:72 residues text "D72" color black height 0.7 font Arial
label #1/A:73 residues text "N73" color black height 0.7 font Arial
label #1/A:74 residues text "S74" color black height 0.7 font Arial
label #1/A:75 residues text "K75" color black height 0.7 font Arial

# Offset labels for readability
label #1/A:71,74,75 residues offset 1,1,0.5
label #1/A:73 residues offset -0.5,1,0.5
label #1/A:72 residues offset -1,0.5,0.5

# Show hydrogen bonds from sites 73-75 to any protein residue (cyan dashed)
hbonds (#1/A:73-75 & ~solvent) restrict (#1 & ~solvent) reveal true color cyan dashes 5

# Show H-bond partner residues as gray ball-and-stick
show #1/A:29,32,52,77 atoms
style #1/A:29,32,52,77 ball
color #1/A:29,32,52,77 gray atoms
color #1/A:29,32,52,77 byhet atoms

# Hide solvent and ions
hide solvent atoms
hide ions

# Center view on the region of interest
view #1/A:71-75
```

## 8G3Z — IGHV1 (FNI17 Fab + Neuraminidase)

Note: PDB residue numbering is offset by +1 from Chothia numbering for this structure. Chothia sites 71-75 correspond to PDB residues 72-76 on chain E.

```
# Load structure
open 8g3z

# Set white background
set bgColor white

# Hide all atom representations, show cartoon ribbons
hide #1 atoms
cartoon #1

# Show only one heavy chain copy (chain E)
hide #1/A,B,C,D,F,G,H,I,J,K,L cartoons

# Color the heavy chain ribbon gray
color #1/E gray cartoons

# Set simple lighting with no depth cueing
lighting simple
lighting depthCue false

# Add silhouette outlines for clarity
graphics silhouettes true width 1

# Increase rendering quality
graphics quality 4

# Show Chothia sites 73-75 (K73, S74, T75 = PDB 74-76) as yellow ball-and-stick
show #1/E:74-76 atoms
style #1/E:74-76 ball
color #1/E:74-76 yellow atoms
color #1/E:74-76 byhet atoms

# Show Chothia site 72 (D72 = PDB 73) as gray ball-and-stick
show #1/E:73 atoms
style #1/E:73 ball
color #1/E:73 gray atoms
color #1/E:73 byhet atoms

# Label residues with single-letter codes (Chothia numbering)
label #1/E:73 residues text "D72" color black height 0.7 font Arial
label #1/E:74 residues text "K73" color black height 0.7 font Arial
label #1/E:75 residues text "S74" color black height 0.7 font Arial
label #1/E:76 residues text "T75" color black height 0.7 font Arial

# Offset labels for readability
label #1/E:75-76 residues offset 1,1,0.5
label #1/E:74 residues offset -0.5,1,0.5
label #1/E:73 residues offset -1,0.5,0.5

# Show hydrogen bonds from sites 73-75 to any residue (cyan dashed)
hbonds #1/E:74-76 restrict any reveal true color cyan dashes 5

# Show H-bond partner residues as gray ball-and-stick
show #1/E:23,78 atoms
style #1/E:23,78 ball
color #1/E:23,78 gray atoms
color #1/E:23,78 byhet atoms

# Hide solvent and ions
hide solvent
hide ions

# Center view on the region of interest
view #1/E:72-76
```
