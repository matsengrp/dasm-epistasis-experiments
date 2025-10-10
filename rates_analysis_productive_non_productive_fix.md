# Mathematical Documentation: Mutation Rate Validation Analysis

**Notebook**: `dasm-epistasis-experiments/rates_analysis_productive_non_productive_fix.ipynb`

## Overview

This notebook implements a count-based validation of DASM and DNSM selection factors by computing empirical mutation rate ratios. The core mathematical principle is that selection factors should equal the ratio of observed (productive) to expected (neutral, non-productive) mutation rates, providing an independent validation of model predictions.

## Mathematical Framework

### 1. Fundamental Selection Factor Relationship

The analysis exploits the following decomposition of selection factors:

$$s = \frac{\mu_{\text{observed}}}{\mu_{\text{expected}}}$$

where:
- $s$ represents a selection factor (multiplicative effect of selection on mutation probability)
- $\mu_{\text{observed}}$ is the mutation rate in productive (functional) BCR sequences under selection
- $\mu_{\text{expected}}$ is the neutral mutation rate from non-productive (out-of-frame) sequences

### 2. Mutation Rate Computation

#### 2.1 Site-Level Mutation Rate (DNSM Comparison)

For a nucleotide site $\ell$ in V-gene family $\mathcal{V}$, the mutation rate is:

$$\mu_{\ell}^{\mathcal{V}} = \frac{N_{\ell}^{\text{mut}}}{\sum_{b \in \mathcal{B}_{\ell}^{\mathcal{V}}} L_b}$$

where:
- $N_{\ell}^{\text{mut}} = \sum_{b \in \mathcal{B}_{\ell}^{\mathcal{V}}} \mathbb{1}[\text{mutation at } \ell \text{ on branch } b]$
- $\mathcal{B}_{\ell}^{\mathcal{V}}$ is the set of phylogenetic parent-child pairs (PCPs) in V-family $\mathcal{V}$ where site $\ell$ has germline identity in the parent sequence
- $L_b$ is the branch length or evolutionary distance for branch $b$

#### 2.2 Branch Length Normalization

Two branch length measures are computed:

**Traditional branch length** (from IQTree phylogenetic inference):
$$L_b^{\text{tree}} = \text{branch\_length}_b$$

**Synonymous mutation frequency** (preferred for consistency with thrifty):
$$L_b^{\text{syn}} = \frac{M_b^{\text{syn}}}{3 \cdot K_b}$$

where:
- $M_b^{\text{syn}} = \sum_{k=1}^{K_b} m_k^{\text{syn}}$ is the total synonymous nucleotide mutations on branch $b$
- $m_k^{\text{syn}}$ is the number of synonymous nucleotide differences in codon $k$
- $K_b$ is the number of codons on branch $b$
- The factor of 3 converts to nucleotide sequence length

The synonymous mutation count for a codon is:

$$m_k^{\text{syn}} = \begin{cases}
d_H(c_{\text{parent}}, c_{\text{child}}) & \text{if } \text{AA}(c_{\text{parent}}) = \text{AA}(c_{\text{child}}) \\
0 & \text{otherwise}
\end{cases}$$

where $d_H$ is the Hamming distance and $\text{AA}(c)$ translates codon $c$ to its amino acid.

#### 2.3 Amino Acid-Specific Mutation Rate (DASM Comparison)

For amino acid transitions $i \to j$ at site $\ell$ in V-family $\mathcal{V}$:

$$\mu_{\ell, i \to j}^{\mathcal{V}} = \frac{N_{\ell, i \to j}}{\sum_{b \in \mathcal{B}_{\ell,i}^{\mathcal{V}}} L_b^{\text{syn}}}$$

where:
- $N_{\ell, i \to j}$ counts observed transitions from parent amino acid $i$ to child amino acid $j$ at site $\ell$
- $\mathcal{B}_{\ell,i}^{\mathcal{V}}$ restricts to branches where:
  - Parent sequence has amino acid $i$ at site $\ell$
  - $i$ matches the germline amino acid for the V-gene at site $\ell$

**Mutation counting**:
$$N_{\ell, i \to j} = \sum_{b \in \mathcal{B}_{\ell,i}^{\mathcal{V}}} \mathbb{1}[\text{child}_b(\ell) = j \land m_{\ell}(b) > 0]$$

where $m_{\ell}(b)$ is the number of nucleotide differences in the codon at site $\ell$ on branch $b$.

#### 2.4 Codon-Specific Mutation Rate

At the finest resolution, for codon transitions $c \to c'$ at site $\ell$:

$$\mu_{\ell, c \to c'}^{\mathcal{V}} = \frac{N_{\ell, c \to c'}}{\sum_{b \in \mathcal{B}_{\ell,c}^{\mathcal{V}}} L_b^{\text{syn}}}$$

with the restriction $d_H(c, c') = 1$ (single nucleotide substitution).

**Codon mutation indicator**:
$$N_{\ell, c \to c'} = \sum_{b \in \mathcal{B}_{\ell,c}^{\mathcal{V}}} \mathbb{1}[\text{child\_codon}_b(\ell) = c' \land m_{\ell}(b) > 0]$$

### 3. Selection Factor Ratio Calculation

#### 3.1 Rate Ratio with Pseudocounts

To handle zero-rate edges, a pseudocount is applied:

$$R_{\text{ratio}} = \frac{\mu^{\text{obs}} + \epsilon}{\mu^{\text{exp}} + \epsilon}$$

where $\epsilon = 10^{-3}$ prevents division by zero and stabilizes ratios for rare mutations.

#### 3.2 Log-Transformation and Clipping

For visualization and regression, log-transformed ratios are used:

$$\log R_{\text{ratio}} = \log\left(\frac{\mu^{\text{obs}} + \epsilon}{\mu^{\text{exp}} + \epsilon}\right)$$

clipped to the interval $[-4, 4]$ to prevent outlier domination in statistical analyses.

Similarly, DASM/DNSM selection factors are log-transformed:

$$\log s_{\text{model}} = \log(s^{\text{DASM/DNSM}})$$

with identical clipping.

### 4. Validation Through Linear Regression

The validation hypothesis is:

$$\log s_{\text{model}} \approx \log R_{\text{ratio}}$$

This is tested via ordinary least squares regression:

$$\log s_{\text{model}} = \beta_0 + \beta_1 \log R_{\text{ratio}} + \varepsilon$$

where $\varepsilon \sim \mathcal{N}(0, \sigma^2)$.

**Perfect agreement** would yield $\beta_0 = 0$ and $\beta_1 = 1$.

**Goodness of fit** is quantified by:

$$R^2 = 1 - \frac{\sum_i (\log s_{\text{model},i} - \hat{y}_i)^2}{\sum_i (\log s_{\text{model},i} - \bar{y})^2}$$

where $\hat{y}_i = \beta_0 + \beta_1 \log R_{\text{ratio},i}$ are fitted values.

## Data Processing Pipeline

### 5. Non-Productive Sequence Processing

#### 5.1 Parent-Child Pair Encoding

For each PCP $b$, mutations are encoded at the nucleotide level:

$$\mathbf{m}_b = (\mathbb{1}[\text{parent}_b[i] \neq \text{child}_b[i]])_{i=1}^{L_b}$$

where $L_b$ is the sequence length for branch $b$.

#### 5.2 Validity Masking

A mask vector indicates valid positions:

$$\mathbf{v}_b = (\mathbb{1}[\text{child}_b[i] \neq \texttt{N}])_{i=1}^{L_b}$$

Only positions with $v_b[i] = 1$ are included in analyses.

#### 5.3 Codon-Level Aggregation

Nucleotide-level data is aggregated to codons. For codon position $k = \lfloor i/3 \rfloor + 1$:

$$c_{\text{parent}}[k] = \text{parent}_b[3k-3 : 3k]$$
$$c_{\text{child}}[k] = \text{child}_b[3k-3 : 3k]$$

#### 5.4 Germline Alignment (Non-Productive Data)

Non-productive sequences have out-of-frame VDJ junctions but in-frame V-gene regions. To enable site-level comparison:

1. **IMGT numbering** is applied based on V-gene alignment
2. **Germline annotations** are merged based on $(v\_gene, site)$ pairs
3. **Site indices** in non-productive data are converted from sequential rank to IMGT position

**Rank to IMGT mapping**:
$$\text{site}_{\text{IMGT}} = f(v\_gene, \text{rank})$$

where $f$ is defined by the germline reference for each V-gene.

### 6. Productive Sequence Processing

Productive data uses DNSM/DASM model outputs (OEPlotter):

#### 6.1 Site Substitution Probabilities DataFrame

Contains per-site information for each PCP:
- `site`: IMGT-numbered position
- `pcp_index`: Phylogenetic branch identifier
- `parent_aa`, `child_aa`: Amino acid states
- `parent_codon`, `child_codon`: Codon states
- `branch_length`: Phylogenetic branch length
- `selection_factor`: DASM/DNSM predicted selection factor
- `prob`: Model-predicted mutation probability
- `neutral_prob`: Neutral mutation probability

#### 6.2 Germline Identity Filtering

Critical filter: Only analyze branches where parent has germline identity:

$$\mathcal{B}_{\ell, \text{analyzed}} = \{b \in \mathcal{B}_{\ell} : \text{parent}_b(\ell) = \text{germline}(\ell, v\_gene_b)\}$$

This ensures:
1. Mutation is "away from germline" (first-hit mutations)
2. Consistency with non-productive data (always germline in V-region)
3. Proper selection factor interpretation

#### 6.3 Nucleotide Mutation Counting

For each codon at site $\ell$ on branch $b$:

$$m_{\ell}(b) = \sum_{i=0}^{2} \mathbb{1}[\text{parent\_codon}_b(\ell)[i] \neq \text{child\_codon}_b(\ell)[i]]$$

This enables distinction between:
- Synonymous mutations: $m_{\ell}(b) > 0$ and $\text{parent\_aa}_b(\ell) = \text{child\_aa}_b(\ell)$
- Non-synonymous mutations: $m_{\ell}(b) > 0$ and $\text{parent\_aa}_b(\ell) \neq \text{child\_aa}_b(\ell)$

### 7. Entrenchment Site Analysis

#### 7.1 Entrenchment Classification

Entrenchment sites are loaded from prior DASM analysis. Two categories:

**Within V-family entrenchment**: Sites $(s, \mathcal{V}, i, j)$ where amino acid $i$ is preferentially retained over $j$ within V-family $\mathcal{V}$ at site $s$.

**Between V-family entrenchment**: Sites $(s, i, j)$ showing differential amino acid preferences across V-families.

#### 7.2 Highlighting in Visualizations

Entrenched sites are overlaid on scatter plots comparing $\log R_{\text{ratio}}$ vs $\log s_{\text{model}}$. These sites are expected to show:

$$\mu^{\text{exp}} \gg \mu^{\text{obs}} \implies R_{\text{ratio}} \ll 1 \implies \log R_{\text{ratio}} < 0$$

indicating strong purifying selection despite high mutability.

## Visualization and Analysis

### 8. Rate Trajectory Plots (Per V-Family)

For each V-family $\mathcal{V} \in \{\text{IGHV1}, \text{IGHV3}, \text{IGHV4}\}$:

**Left panel** shows $\mu_{\ell}^{\text{obs}}$ and $\mu_{\ell}^{\text{exp}}$ as functions of IMGT site $\ell$:

$$\begin{aligned}
\mu_{\ell}^{\text{obs}}(\mathcal{V}) &= \frac{N_{\ell}^{\text{mut, productive}}}{\sum_{b \in \mathcal{B}_{\ell}^{\mathcal{V}, \text{productive}}} L_b^{\text{syn}}} \\
\mu_{\ell}^{\text{exp}}(\mathcal{V}) &= \frac{N_{\ell}^{\text{mut, non-productive}}}{\sum_{b \in \mathcal{B}_{\ell}^{\mathcal{V}, \text{non-productive}}} L_b^{\text{syn}}}
\end{aligned}$$

**Right panel** shows the ratio:

$$R_{\text{ratio}}(\ell, \mathcal{V}) = \frac{\mu_{\ell}^{\text{obs}}(\mathcal{V})}{\mu_{\ell}^{\text{exp}}(\mathcal{V})}$$

**CDR regions** are highlighted with red shading (IMGT positions 27-38, 56-65, 105-110).

### 9. Scatter Plots with Regression

#### 9.1 DNSM Validation (Site-Level)

$$\text{scatter}(\log R_{\text{ratio}}(\ell, \mathcal{V}), \log s_{\ell}^{\text{DNSM}}(\mathcal{V}))$$

for all $(s, \mathcal{V})$ pairs. Entrenched sites highlighted by color (site) and marker shape (V-family).

#### 9.2 DASM Validation (Amino Acid-Level)

$$\text{scatter}(\log R_{\text{ratio}}(\ell, i \to j, \mathcal{V}), \log s_{\ell, i \to j}^{\text{DASM}}(\mathcal{V}))$$

**Filtering**: Only transitions with $N^{\text{obs}} + N^{\text{exp}} > 20$ and single-nucleotide accessibility.

**Single-nucleotide accessibility** for amino acids $i \to j$:

$$\exists c, c' : \text{AA}(c) = i, \text{AA}(c') = j, d_H(c, c') = 1$$

This ensures the transition is mutationally accessible and has sufficient statistical power.

### 10. Synonymous vs. Non-Synonymous Analysis

Using codon-level rates, the log-ratio distributions are compared:

**Synonymous transitions**: $\text{AA}(c) = \text{AA}(c')$

$$\log R_{\text{syn}} = \log\left(\frac{\mu_{c \to c'}^{\text{obs}}}{\mu_{c \to c'}^{\text{exp}}}\right)$$

Expected to center near 0 (neutral evolution).

**Non-synonymous transitions**: $\text{AA}(c) \neq \text{AA}(c')$

$$\log R_{\text{nonsyn}} = \log\left(\frac{\mu_{c \to c'}^{\text{obs}}}{\mu_{c \to c'}^{\text{exp}}}\right)$$

Expected to show broader distribution with negative skew (purifying selection dominates).

**Statistical test** (implicit in visualization):

$$H_0: \mathbb{E}[\log R_{\text{syn}}] = 0 \quad \text{vs.} \quad H_0: \mathbb{E}[\log R_{\text{nonsyn}}] < 0$$

Histograms with kernel density estimation (KDE) overlays visualize these distributions.

## Key Implementation Details

### 11. Germline Identity Filtering

For both observed and expected rates, only branches with germline identity in the parent are analyzed:

$$\text{is\_germline\_codon}_b(\ell) = \mathbb{1}[\text{parent\_codon}_b(\ell) = \text{germline\_codon}(v\_gene_b, \ell)]$$

$$\text{is\_germline\_aa}_b(\ell) = \mathbb{1}[\text{parent\_aa}_b(\ell) = \text{germline\_aa}(v\_gene_b, \ell)]$$

Analyses use:
- Codon-level: Filter by `is_germline_codon == True`
- Amino acid-level: Filter by `is_germline_aa == True`

### 12. Leaf Node Removal

To avoid biases from terminal branches (which may have elevated mutation rates due to sequencing errors or true terminal variation):

$$\mathcal{B}_{\text{analyzed}} = \{b \in \mathcal{B} : \neg \text{child\_is\_leaf}_b\}$$

This removes branches where the child node is a leaf in the phylogenetic tree.

### 13. V-Gene Family Restriction

Primary analyses focus on three major V-families:

$$\mathcal{V}_{\text{analyzed}} = \{\text{IGHV1}, \text{IGHV3}, \text{IGHV4}\}$$

These families have sufficient sample sizes and represent the majority of human BCR repertoire.

## Mathematical Properties and Assumptions

### 14. Assumptions

1. **Non-productive neutrality**: Mutations in non-productive sequences are not subject to selection (out-of-frame VDJ prevents functional BCR expression)

2. **V-region consistency**: The V-gene region in non-productive sequences maintains in-frame alignment with germline (frameshifts occur only at VDJ junction)

3. **Mutational process homogeneity**: The underlying mutational process (SHM machinery) is identical in productive and non-productive sequences

4. **Independence**: Mutations at different sites are independent (known simplification; epistasis violates this)

5. **First-hit approximation**: By filtering to germline identity in parent, we approximate first-hit mutation rates

### 15. Statistical Robustness

**Minimum count threshold**: $N_{\text{total}} = N^{\text{obs}} + N^{\text{exp}} > 20$

This ensures:
- Reasonable binomial confidence intervals: $\text{SE}[\hat{\mu}] \approx \sqrt{\mu(1-\mu)/N}$
- Rate ratio stability: $\text{CV}[R_{\text{ratio}}] \propto 1/\sqrt{N_{\text{total}}}$

**Pseudocount justification**: $\epsilon = 10^{-3}$ is chosen to be:
- Small relative to observed rates (typically $10^{-2}$ to $10^{0}$)
- Large enough to stabilize ratios when one rate is zero

### 16. Interpretation of Results

**Strong correlation** ($R^2 > 0.7$): DASM/DNSM models accurately predict selection factors from neutral mutation patterns

**Slope near 1**: Rate ratios quantitatively match selection factor magnitudes

**Intercept near 0**: No systematic bias between methods

**Entrenched sites on regression line**: Previously identified sites are not outliers, validating prior findings through independent methodology

## Biological Interpretation

### 17. Entrenchment Mechanism

Sites showing high $\mu^{\text{exp}}$ but low $\mu^{\text{obs}}$ indicate:

$$\text{High intrinsic mutability} + \text{Strong purifying selection} \implies \text{Functional constraint}$$

The mathematical signature:

$$\mu^{\text{exp}} \gg \mu^{\text{obs}} \iff s \ll 1 \iff \log s < 0$$

identifies sites under strong negative selection despite being mutational hotspots.

### 18. V-Family Specific Selection

Variation in $R_{\text{ratio}}(\ell, \mathcal{V})$ across V-families indicates differential selection pressures:

$$\text{Var}_{\mathcal{V}}[R_{\text{ratio}}(\ell, \mathcal{V})] > 0 \implies \text{V-family context-dependent selection}$$

This suggests epistatic interactions between site $\ell$ and V-family-specific structural contexts.

### 19. CDR Enrichment

Elevated $\mu^{\text{obs}}$ in CDR regions (despite selection) reflects:

$$\mu^{\text{CDR, obs}} > \mu^{\text{FWR, obs}} \quad \text{due to} \quad \mu^{\text{CDR, exp}} \gg \mu^{\text{FWR, exp}}$$

indicating that even with selection against many mutations, the intrinsically high mutation rate in CDRs leads to observable diversity.

## Conclusion

This notebook provides a rigorous, count-based validation of DASM and DNSM selection factors through empirical mutation rate ratios. The strong agreement between $\log R_{\text{ratio}}$ and model predictions confirms that:

1. DASM/DNSM accurately decompose observed mutation patterns into neutral and selective components
2. Non-productive sequences provide a valid neutral baseline
3. Previously identified entrenched sites show expected mathematical signatures (high neutral rates, strong selection)
4. The methodology is robust across granularities (site, amino acid, codon level)

The mathematical framework establishes a principled connection between population-level mutation frequencies and selection factors, enabling independent validation of machine learning-derived selection landscapes.
