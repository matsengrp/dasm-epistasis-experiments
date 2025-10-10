# Mathematical Analysis Overview: Validation of Entrenchment Sites Using Mutation Rate Ratios

## Summary

This analysis validates entrenchment findings in B-cell receptor (BCR) sequences by comparing observed mutation rates (from productive sequences) to expected mutation rates (from non-productive sequences). The core mathematical framework establishes that the ratio of these rates should equal selection factors predicted by DASM (Detailed Amino acid Selection Model) and DNSM (Detailed Nucleotide Selection Model).

**Key File**: `dasm-epistasis-experiments/rates_analysis_productive_non_productive_fix.ipynb`

## Mathematical Framework

### Selection Factor Decomposition

The fundamental relationship exploited in this analysis is:

$$s_{ij}^{\text{DASM}} = \frac{\mu_{ij}^{\text{obs}}}{\mu_{ij}^{\text{neutral}}}$$

where:
- $s_{ij}^{\text{DASM}}$ is the DASM-predicted selection factor for mutation from amino acid $i$ to $j$
- $\mu_{ij}^{\text{obs}}$ is the observed mutation rate (productive sequences under selection)
- $\mu_{ij}^{\text{neutral}}$ is the expected neutral mutation rate (non-productive sequences)

This relationship extends to the DNSM model at the nucleotide level:

$$s_{\ell}^{\text{DNSM}} = \frac{\mu_{\ell}^{\text{obs}}}{\mu_{\ell}^{\text{neutral}}}$$

where $\ell$ indexes nucleotide sites.

### Mutation Rate Estimation

#### Site-Level Mutation Rate

For a given site $\ell$ in V-family $\mathcal{V}$, the mutation rate is defined as:

$$\mu_{\ell}^{\mathcal{V}} = \frac{N_{\ell}^{\text{mut}}}{\sum_{b \in \mathcal{B}_{\ell}^{\mathcal{V}}} L_b^{\text{syn}}}$$

where:
- $N_{\ell}^{\text{mut}}$ is the number of mutations observed at site $\ell$
- $\mathcal{B}_{\ell}^{\mathcal{V}}$ is the set of phylogenetic branches containing site $\ell$ with germline identity in V-family $\mathcal{V}$
- $L_b^{\text{syn}}$ is the synonymous mutation frequency for branch $b$, serving as a normalized measure of evolutionary time

#### Synonymous Mutation Frequency as Branch Length

The key innovation for rate normalization is using synonymous mutation frequency instead of traditional branch lengths:

$$L_b^{\text{syn}} = \frac{\sum_{k=1}^{K_b} \mathbb{1}[\text{syn}(k)]}{3K_b}$$

where:
- $K_b$ is the number of codons in branch $b$
- $\mathbb{1}[\text{syn}(k)]$ indicates whether codon $k$ contains a synonymous mutation
- The factor of 3 normalizes to nucleotide length

This measure provides consistency with thrifty calculations and controls for variable mutational loads across branches.

#### Amino Acid-Specific Mutation Rate

For DASM validation, rates are stratified by parent and child amino acids. For a transition $i \to j$ at site $\ell$ in V-family $\mathcal{V}$:

$$\mu_{\ell,i \to j}^{\mathcal{V}} = \frac{N_{\ell,i \to j}}{\sum_{b \in \mathcal{B}_{\ell,i}^{\mathcal{V}}} L_b^{\text{syn}}}$$

where:
- $N_{\ell,i \to j}$ is the number of observed transitions from amino acid $i$ to $j$ at site $\ell$
- $\mathcal{B}_{\ell,i}^{\mathcal{V}}$ restricts to branches where site $\ell$ has parent amino acid $i$ matching the germline

#### Codon-Level Mutation Rate

At the finest resolution, codon-specific rates are computed:

$$\mu_{\ell,c \to c'}^{\mathcal{V}} = \frac{N_{\ell,c \to c'}}{\sum_{b \in \mathcal{B}_{\ell,c}^{\mathcal{V}}} L_b^{\text{syn}}}$$

where transitions are restricted to single-nucleotide substitutions: $d_H(c, c') = 1$ (Hamming distance).

### Validation Methodology

#### Log-Ratio Comparison

To validate DASM/DNSM predictions, we analyze the log-transformed rate ratio:

$$\log s_{\text{ratio}} = \log\left(\frac{\mu^{\text{obs}} + \epsilon}{\mu^{\text{neutral}} + \epsilon}\right)$$

where $\epsilon = 10^{-3}$ is a pseudocount to handle zero rates. This is compared against model predictions:

$$\log s_{\text{model}} = \log(s^{\text{DASM/DNSM}})$$

Values are clipped to $[-4, 4]$ to prevent outlier domination.

#### Linear Regression Analysis

The validation quantifies agreement through linear regression:

$$\log s_{\text{model}} = \beta_0 + \beta_1 \log s_{\text{ratio}} + \varepsilon$$

where agreement is assessed via:
- Slope $\beta_1$ (should approximate 1 under perfect agreement)
- $R^2$ coefficient (measures explained variance)

### Entrenchment Site Highlighting

The analysis specifically highlights previously identified entrenched sites, characterized by:

1. **Between-family entrenchment**: Sites showing differential selection across V-families
2. **Within-family entrenchment**: Sites showing entrenchment within a specific V-family (more surprising finding)

These sites are expected to show:
- High neutral mutation rates $\mu^{\text{neutral}}$ (mutational hotspots)
- Low observed mutation rates $\mu^{\text{obs}}$ (strong purifying selection)
- Consequently, $s < 1$ (selection against mutation)

### CDR Region Analysis

The analysis incorporates IMGT-numbered CDR regions:
- CDR1: sites 27-38
- CDR2: sites 56-65
- CDR3: sites 105-110 (V-gene portion only)

These regions are highlighted in visualizations as they typically show elevated mutation and selection patterns.

## Statistical Considerations

### Filtering Criteria

To ensure statistical robustness, analyses apply:

1. **Minimum mutation count**: $N^{\text{obs}} + N^{\text{neutral}} > 20$ for amino acid-level analysis
2. **Single-mutation accessibility**: Only amino acid pairs reachable by single nucleotide substitution
3. **Germline restriction**: Only branches with germline identity at the analyzed site
4. **Non-productive sequence validity**: Only V-gene encoded positions (VDJ junction may be out-of-frame)

### Synonymous vs. Non-synonymous Rates

The codon-level analysis enables decomposition:

$$\mu_{\ell,c \to c'}^{\mathcal{V}} = \begin{cases}
\mu_{\ell}^{\text{syn}} & \text{if } \text{AA}(c) = \text{AA}(c') \\
\mu_{\ell}^{\text{nonsyn}} & \text{if } \text{AA}(c) \neq \text{AA}(c')
\end{cases}$$

Synonymous mutations serve as a neutral baseline, while non-synonymous rates reflect selection.

## Key Findings

1. **DNSM validation**: Strong correlation between log-ratio and DNSM selection factors ($R^2 \approx$ reported in plots)
2. **DASM validation**: Robust agreement for amino acid-specific transitions after filtering
3. **Entrenchment consistency**: Previously identified entrenched sites (38, 40, 55, 57, 66) show:
   - High $\mu^{\text{neutral}}$ (confirming mutational hotspot status)
   - Strong depletion in $\mu^{\text{obs}}$ (confirming selection)
   - Consistent with DASM/DNSM predictions (not outliers)

4. **Synonymous vs. non-synonymous**: Clear separation in log-ratio distributions, with non-synonymous showing broader distribution due to selection

## Implementation Notes

- **Data sources**:
  - Observed rates from productive sequences (DNSM/DASM OEPlotter output)
  - Expected rates from non-productive out-of-frame sequences (tangshm dataset)
- **Branch filtering**: Leaf nodes removed to avoid terminal branch artifacts
- **IMGT alignment**: Non-productive sequences aligned to germline V-gene for consistent site numbering
