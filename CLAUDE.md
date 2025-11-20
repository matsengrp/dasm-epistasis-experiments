# Project Context for Claude

## Project Overview

This repository (`dasm-epistasis-experiments`) focuses on modeling antibody evolution using deep learning approaches that separate mutation and selection processes:

**Core Concept**: Antibodies evolve through somatic hypermutation (SHM) and natural selection in germinal centers. Traditional language models conflate these processes, but we model them separately:

- **Mutation Model**: Captures the neutral process of somatic hypermutation using nucleotide-level biases
- **Selection Model**: Uses transformer networks to predict which amino acid changes are favored or disfavored by natural selection

**Two Main Approaches**:
1. **DNSM (Deep Natural Selection Model)**: Predicts a single selection factor per site indicating whether that site experiences purifying (<1) or diversifying (>1) selection (output_dim=1)
2. **DASM (Deep Amino Acid Selection Model)**: Predicts selection factors for every possible amino acid substitution at every site, enabling direct comparison with functional assays (output_dim=20)

Both models use the same underlying neural network architecture but are trained differently on parent-child pairs from reconstructed phylogenetic trees of antibody clonal families.

## Environment Setup

To activate the correct Python environment:

```bash && 
source ~/.bashrc && conda activate netam_env
```

## Development Principles

### Fail-Fast, No Fallbacks
- **No Silent Fallbacks**: Code must fail immediately when expected conditions aren't met. Silent fallback behavior masks bugs and creates unpredictable systems.
- **Explicit Error Messages**: When something goes wrong, stop execution with clear error messages explaining what failed and what was expected.
- **Example**: `raise ValueError(f"Required model {model_name} not found")` instead of falling back to first available model.

## Development Tools

### GitHub Integration

Use `gh` CLI for issue management:

```bash
gh issue create --title "Title" --body-file issue.md
gh issue list
gh pr create --title "Title" --body "Description"
```

## Important Note on File Locations

**Papers**: The `papers/` directory contains reference PDFs related to the research. These can be read using MCP PDF tools.

## Planning Documents

When creating planning documents for complex tasks or analysis workflows, use the `*PLAN.md` suffix (e.g., `ANALYSIS_PLAN.md`, `MODEL_SCORING_PLAN.md`). These files are automatically gitignored as working documents and should contain:

- **Task breakdown**: Step-by-step implementation plans with code snippets
- **Status tracking**: What's completed, in progress, or pending
- **Infrastructure decisions**: Architecture choices and integration strategies
- **Expected outputs**: Clear deliverables and success criteria

This keeps planning work organized while avoiding clutter in version control.
