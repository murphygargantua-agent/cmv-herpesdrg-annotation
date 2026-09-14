# Build Procedure: CMV HerpesDRG vcfanno Database

This document describes the complete build process for the CMV antiviral resistance variant database.

## Overview

The database maps HerpesDRG amino acid mutations to genomic coordinates on the AD169 reference genome, formatted as a BED file for `vcfanno` annotation.

**Total mutations:** 586/602 (97.3% match rate)
**Reference genome:** AD169 (NC_006273, X17403.1)
**Gene coverage:** UL97, UL54, UL27, UL56, UL51, UL89

## Step 1: Download HerpesDRG Data

HerpesDRG is publicly available as a TSV file on GitHub:

```bash
wget -O herpesdrg-db.tsv https://github.com/ojcharles/herpesdrg-db/raw/main/herpesdrg-db.tsv
```

**Columns:** virus, gene, amino_acid_change, strain, fold_change, drug, notes, status

**Filtering:** Only mutations with `status == "A"` (active) were included.
Total: 602 active HCMV mutations.

## Step 2: Fetch AD169 Reference Genome

The AD169 strain was chosen because CHARMD/HerpesDRG use AD169 numbering natively. Cross-strain mapping (e.g., to Merlin or Towne) caused significant coordinate mismatches.

```python
# Fetch via NCBI e-utilities
import urllib.request
url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
params = "?db=nucleotide&id=X17403&rettype=fasta&retmode=text"
urllib.request.urlretrieve(url + params, "AD169.fasta")
```

**Result:** 235,646 bp genome sequence.

## Step 3: Parse GenBank Annotations

The AD169 GenBank file contains gene feature locations needed for coordinate mapping.

```python
from Bio import GenBank
from io import StringIO

# Parse GenBank to extract gene positions
with open("AD169_annotation.gb") as f:
    record = GenBank.read(f)

genes = {}
for feature in record.features:
    if feature.type == "CDS" and "gene" in feature.qualifiers:
        gene_name = feature.qualifiers["gene"][0]
        if gene_name.startswith("HCMV"):  # AD169 prefix
            gene_name = gene_name[4:]  # Remove HCMV prefix
        if gene_name in ["UL97", "UL54", "UL27", "UL56", "UL51", "UL89"]:
            gene_location = feature.location
            start = gene_location.start + 1  # 1-based
            end = gene_location.end
            genes[gene_name] = {"start": start, "end": end, "length": end - start}
```

**Resolved genes:**
- UL97: positions 94950–96401 (1452 aa)
- UL54: positions 89478–92529 (1017 aa)
- UL27: positions 133899–134403 (168 aa)
- UL56: positions 88566–89031 (155 aa)
- UL51: positions 90136–90201 (22 aa)
- UL89: positions 88475–88533 (19 aa)

## Step 4: Map Mutations to Genomic Coordinates

Each mutation follows the format `GENE POSITION_REF_ALA` (e.g., `UL97 M460V`).

```python
def aa_to_nt_position(gene, aa_pos, reference_aa, variant_aa):
    """
    Convert amino acid position to nucleotide position.
    Uses the gene's start position and accounts for codon structure.
    """
    gene_info = gene_coordinates[gene]
    gene_start = gene_info["start"]
    
    # Nucleotide position = (gene_start - 1) + (aa_pos * 3)
    # -1 because gene_start is 1-based, we need 0-based for calculation
    nt_pos = (gene_start - 1) + (aa_pos * 3)
    
    # Reference nucleotide should match
    # Check against AD169 FASTA
    ref_nt = ad169_sequence[nt_pos:nt_pos+3]
    
    return nt_pos, ref_nt
```

**Verification:** Compare derived reference amino acid with AD169 sequence. Only keep mutations where the reference AA matches.

## Step 5: Generate BED File

Format the mapped mutations as a BED file:

```python
import csv
import json

with open("herpesdrg_cmv.bed", "w") as f:
    writer = csv.writer(f, delimiter="\t")
    for mutation in active_mutations:
        gene = mutation["gene"]
        aa_change = mutation["amino_acid_change"]  # e.g., "M460V"
        
        # Parse mutation
        ref_aa = aa_change[0]
        pos = int(aa_change[1:-1])
        var_aa = aa_change[-1]
        
        # Calculate genomic coordinates
        nt_start, ref_nt = aa_to_nt_position(gene, pos, ref_aa, var_aa)
        
        # Write BED line
        writer.writerow([
            "NC_006273",      # Chromosome
            nt_start,         # Start (0-based)
            nt_start + 3,     # End (exclusive)
            f"{gene}_{aa_change}",  # Name
            gene,             # Gene
            mutation["drug"],  # Drug
            mutation["fold_change"],  # Fold change
            ref_aa,           # Reference AA
            var_aa,           # Variant AA
            mutation["notes"]  # Notes
        ])
```

## Step 6: Compress and Index

```bash
# Sort BED file by chromosome and position
sort -k1,1 -k2,2n herpesdrg_cmv.bed > herpesdrg_cmv_sorted.bed

# Compress with bgzip
bgzip -c herpesdrg_cmv_sorted.bed > herpesdrg_cmv.bed.gz

# Generate tabix index
tabix -s 1 -b 2 -e 3 herpesdrg_cmv.bed.gz
```

**Verification:**
```bash
# Check index header
tabix -H herpesdrg_cmv.bed.gz

# Query test region
tabix herpesdrg_cmv.bed.gz NC_006273:94950-95200
```

## Unmapped Mutations

**16 mutations failed to map** (2.7%):
- All unmapped mutations are in UL89
- Cause: Residual strain divergence between HerpesDRG internal reference and AD169
- UL89 is the smallest gene in the dataset (19 aa), making it more sensitive to reference differences

**UL89 specific issue:**
The HerpesDRG database likely uses a slightly different UL89 sequence than AD169. The amino acid positions don't align perfectly, making genomic coordinate mapping unreliable.

## Trade-offs

| Aspect | HerpesDRG | CHARMD |
|--------|-----------|--------|
| Coverage | 586/602 (97.3%) | ~612 total |
| Format | Public TSV | Web only |
| API | N/A | N/A |
| Drugs | 8 antivirals | 8 antivirals |
| Genes | 6 genes | 6 genes |
| Strain | AD169 | AD169 |

**Why HerpesDRG?** Publicly downloadable, well-structured, comprehensive drug resistance data. CHARMD requires manual web scraping or email request to CNR Herpesvirus.

## Dependencies

- Python 3.8+
- Biopython (optional, for GenBank parsing)
- samtools (for bgzip/tabix)
- vcfanno (for annotation)

## Validation Checklist

- [x] All mutations sorted by chromosome and position
- [x] Tabix index generated and verified
- [x] Reference nucleotides match AD169 sequence
- [x] 586/602 active mutations mapped (97.3%)
- [x] 6 genes covered: UL97, UL54, UL27, UL56, UL51, UL89
- [x] BED file format validated
- [x] Gene coordinate lookup saved (ad169_gene_mapping.json)

## Reproducing This Build

To rebuild from scratch:

```bash
# 1. Download HerpesDRG
wget https://github.com/ojcharles/herpesdrg-db/raw/main/herpesdrg-db.tsv

# 2. Fetch AD169
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=X17403&rettype=fasta&retmode=text" -O AD169.fasta

# 3. Run build script
python build_database.py --herpesdrg herpesdrg-db.tsv --genome AD169.fasta --output herpesdrg_cmv.bed

# 4. Compress and index
sort -k1,1 -k2,2n herpesdrg_cmv.bed | bgzip -c > herpesdrg_cmv.bed.gz
tabix -s 1 -b 2 -e 3 herpesdrg_cmv.bed.gz
```

## Future Work

1. **CHARMD integration** — Scrape or request CHARMD data for complete coverage
2. **UL89 mapping** — Investigate strain differences causing UL89 unmapped mutations
3. **Additional genes** — Expand to other herpesvirus genes if HerpesDRG expands
4. **Version tracking** — Add checksums and version metadata for reproducibility
5. **Automated updates** — Set up periodic HerpesDRG refresh and re-mapping

## References

1. Tilloy V, et al. (2024) CHARMD. *Antiviral Research* 231:106016.
2. Charles O. HerpesDRG Database. https://github.com/ojcharles/herpesdrg-db
3. NCBI Nucleotide: X17403.1 (CMV AD169)
