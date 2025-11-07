# Support Chothia numbering scheme in addition to IMGT

## Background

Currently, the codebase assumes IMGT numbering throughout:
- CDR boundary definitions use IMGT ranges (CDR1: 27-38, CDR2: 56-65)
- ANARCI output processing includes IMGT-specific validation checks
- Multiple components depend on IMGT numbering scheme

## Components that need modification

### 1. CDR boundary definitions

**Location:** `pcp-pipeline/germline/scripts/write_ogrdb_gene_nt_cdrs.py`

Currently downloads CDR1 and CDR2 boundaries from OGRDB using IMGT numbering (lines 24-27):
```python
imgt_cdr1_aa_start = 27
imgt_cdr1_aa_end   = 38
imgt_cdr2_aa_start = 56
imgt_cdr2_aa_end   = 65
```

CDR3 boundaries are extracted from partis annotations (V(D)J junction).

**Change needed:** Create parallel Chothia version that downloads/defines Chothia CDR1/CDR2 boundaries from OGRDB. CDR3 extraction from partis can remain the same.

### 2. ANARCI file processing in `netam/netam/oe_plot.py`

**Location:** `netam/netam/oe_plot.py` - `get_numbering_dict()` function (lines 1190-1296)

**Current issues:**
- IMGT-specific validation checks (lines 1230-1244):
  - Rejects insertions other than `111.*` or `112.*`
  - Validates CDR annotation against IMGT numbering

**Changes needed:**
- Start by removing all validation checks (can add Chothia-specific checks later if needed)
- Make numbering scheme configurable

### 3. DASM/DNSM `write_sites_oe` function

**Location:** DASM/DNSM oe_plot utilities (e.g., `dxsm_oe.py`)

**Current issue:**
- `write_sites_oe` looks up ANARCI files in a dictionary that only contains IMGT-aligned files
- No option to specify a Chothia ANARCI file

**Changes needed:**
- Add an additional Chothia ANARCI file dictionary to DASM/DNSM data containers
- Add parameter to `write_sites_oe` to specify which numbering scheme dictionary to use
- Detect numbering scheme from filename and set accordingly (e.g., if "chothia" in filename)

### 4. CDR annotation functions

**Location:** `netam/netam/oe_plot.py`

Functions affected:
- `is_imgt_cdr()` (lines 1164-1187) - hardcoded IMGT ranges
- `pcp_sites_cdr_annotation()` (lines 1299-1330) - relies on PCP CDR boundaries

**Changes needed:**
- Generalize `is_imgt_cdr()` → `is_cdr(site, scheme='imgt')`
- Make scheme selection configurable

### 5. PCP file generation

**Location:** `pcp-pipeline/python/utils.py`

Function `get_cdr_codon_positions()` (lines 436-470) looks up CDR boundaries from CSV files.

**Changes needed:**
- Add `numbering_scheme` parameter to `write_pcp_df()` and related functions
- Select appropriate CDR CSV file based on scheme parameter (CDR1/CDR2 only; CDR3 unchanged)

### 6. ANARCI pipeline

**Changes needed:**
- Run ANARCI with `--scheme chothia` option to generate Chothia-numbered output files
- Store both IMGT and Chothia ANARCI output files separately
- Data containers need dictionaries for both IMGT and Chothia ANARCI files

## Proposed implementation

1. **Add numbering scheme parameter** throughout the pipeline:
   ```python
   numbering_scheme='imgt'  # or 'chothia'
   ```

2. **Create Chothia CDR definition files:**
   - New script or modified script: Download/define Chothia CDR1/CDR2 boundaries from OGRDB
   - Generate: `ogrdb_*_cdr1_cdr2_chothia.csv` files
   - Note: Chothia boundaries are more variable and length-dependent than IMGT
   - CDR3 definitions remain unchanged (extracted from partis)

3. **Make ANARCI processing scheme-aware:**
   - In `get_numbering_dict()` (netam/oe_plot.py): Remove IMGT validation checks initially
   - Can add Chothia-specific validation later if needed
   - Generalize CDR checking: `is_cdr(site, scheme='imgt')`

4. **Update DASM/DNSM data containers:**
   - Add Chothia ANARCI file dictionaries alongside existing IMGT dictionaries
   - Modify `write_sites_oe` to accept numbering scheme parameter
   - Auto-detect scheme from filename if "chothia" appears in path

5. **Run ANARCI pipeline for Chothia:**
   - Execute ANARCI with `--scheme chothia` for all datasets
   - Store Chothia output files separately from IMGT files
   - Maintain both sets of files for flexibility

## Files requiring changes

- `pcp-pipeline/germline/scripts/write_ogrdb_gene_nt_cdrs.py` - create Chothia version for CDR1/CDR2
- `pcp-pipeline/python/utils.py` - add scheme parameter to PCP functions
- `netam/netam/oe_plot.py`:
  - `get_numbering_dict()` - remove IMGT-specific validation checks
  - `is_imgt_cdr()` - generalize to support multiple schemes
- DASM oe_plot utilities (e.g., `dxsm_oe.py`):
  - Add Chothia ANARCI file dictionary to data containers
  - Add numbering scheme parameter to `write_sites_oe`
- DNSM oe_plot utilities:
  - Add Chothia ANARCI file dictionary to data containers
  - Add numbering scheme parameter to `write_sites_oe`
- ANARCI pipeline scripts - run with `--scheme chothia` option

## Additional considerations

- **Chothia complexity:** Chothia CDR boundaries are more variable and length-dependent than IMGT, may need V gene family-specific logic
- **Backward compatibility:** Default to IMGT to maintain existing behavior
- **ANARCI runs:** Need to run ANARCI with `--scheme chothia` for all datasets to generate Chothia-numbered files
- **Validation:** Different numbering schemes may require different validation rules

## References

- IMGT numbering: https://www.imgt.org/IMGTScientificChart/Numbering/IMGTnumberingCDR_VH.html
- Chothia numbering: Chothia & Lesk (1987), Al-Lazikani et al. (1997)
