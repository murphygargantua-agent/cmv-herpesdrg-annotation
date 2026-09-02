# CMV Antiviral Resistance Variant Database

A curated database of cytomegalovirus (CMV) antiviral resistance mutations formatted for integration with `vcfanno`. Built from [HerpesDRG](https://github.com/ojcharles/herpesdrg-db) data, mapped to the [AD169](https://www.ncbi.nlm.nih.gov/nuccore/X17403) reference genome.

## Quick Start

### Using with vcfanno

```bash
# Install vcfanno
# https://github.com/brentp/vcfanno

vcfanno vcfanno_cmv.toml your_variants.vcf.gz > annotated_variants.vcf
```

**vcfanno_cmv.toml:**
```toml
[[annotation]]
file = "data/herpesdrg_cmv.bed.gz"
fields = ["name"]
names = ["cmv_resistance"]
ops = ["self"]
```

### Output Format

Annotations appear in the VCF INFO field as:
```
UL97_M460V|Ganciclovir=5.2;Cidofovir=1.8
UL54_D542E|Foscarnet=23.0;Cidofovir=12.0
```

## Database Contents

| Metric | Value |
|--------|-------|
| Total mutations | 586 (97.3% match rate) |
| Genes covered | UL97, UL54, UL27, UL56, UL51, UL89 |
| Reference genome | AD169 (NC_006273) |
| Source | HerpesDRG (LANE-CMR) |
| Drugs covered | 8 antivirals |

### Gene Distribution

| Gene | Count | Drug Classes |
|------|-------|--------------|
| UL54 | 270 | Ganciclovir, Foscarnet, Cidofovir, Brincidofovir, Valganciclovir |
| UL97 | 184 | Ganciclovir, Valganciclovir, Maribavir |
| UL56 | 104 | Letermovir |
| UL51 | 9 | Letermovir |
| UL27 | 19 | Maribavir |
| UL89 | — | Letermovir |

## Data Sources

### Primary: HerpesDRG (LANE-CMR)

- **Source:** [HerpesDRG Database](https://github.com/ojcharles/herpesdrg-db)
- **Download:** `herpesdrg-db.tsv`
- **Coverage:** 602 HCMV mutations (586 mapped = 97.3%)
- **Format:** Tab-separated with gene, amino acid change, fold-change values, drug associations

### Reference: AD169 Genome

- **GenBank:** [NC_006273](https://www.ncbi.nlm.nih.gov/nuccore/X17403)
- **NCBI Accession:** X17403.1
- **Length:** 235,646 bp
- **Why AD169?** CHARMD/HerpesDRG use AD169 strain numbering natively. Cross-strain mapping to Merlin or Towne caused significant coordinate mismatches.

## Build Process

### Step-by-Step

1. **Download HerpesDRG TSV** from GitHub
2. **Fetch AD169 genome** via NCBI e-utilities (`efetch.fcgi`)
3. **Parse GenBank features** to extract gene coordinates
4. **Map amino acid changes** to genomic coordinates:
   - Parse mutation format (e.g., `UL97 M460V`)
   - Calculate nucleotide position from amino acid position
   - Verify reference nucleotide matches AD169 sequence
5. **Generate BED file** with chrom, start, end, and metadata
6. **Compress and index** with `bgzip` and `tabix`

### Known Limitations

- **16 unmapped mutations** (2.7%) — primarily UL89 due to residual strain divergence
- **UL89 coverage** is incomplete compared to other genes
- **Fold-change data** may have missing values (represented as `N/A`)
- **Drug associations** are based on HerpesDRG annotations which may lag behind recent publications

### CHARMD Note

The [CHARMD](https://www.unilim.fr/cnr-herpesvirus/outils/codexmv/) database (Tilloy et al., 2024) is the gold standard but lacks programmatic access (no API, no bulk download). HerpesDRG was chosen for its public availability and comprehensive coverage. CHARMD data (~612 mutations) is available via manual web scraping but has not been integrated here.

## File Structure

```
cmv-herpesdrg-annotation/
├── README.md                    # This file
├── BUILD_PROCEDURE.md           # Detailed build methodology
├── vcfanno_cmv.toml             # vcfanno configuration
├── data/
│   ├── herpesdrg_cmv.bed        # Plain text BED format
│   ├── herpesdrg_cmv.bed.gz     # Compressed BED (bgzip)
│   ├── herpesdrg_cmv.bed.gz.tbi # Tabix index
│   └── ad169_gene_mapping.json  # Gene coordinate lookup table
└── scripts/
    ├── build_database.py        # Main build script
    └── fetch_genome.py          # Genome fetching script
```

## BED File Format

Each line contains tab-separated columns:

| Column | Description | Example |
|--------|-------------|---------|
| 1 | Chromosome | `NC_006273` |
| 2 | Start position (0-based) | `94961` |
| 3 | End position | `94964` |
| 4 | Mutation name | `UL97_M460V` |
| 5-12 | Gene, drug, fold-changes, notes | ... |

### Columns

1. **chrom** — NC_006273 (AD169)
2. **start** — 0-based start position
3. **end** — 1-based end position (exclusive)
4. **name** — Gene_AminoAcidChange (e.g., `UL97_M460V`)
5. **gene** — Gene name
6. **drug** — Associated drug
7. **fold_change** — EC50 fold change value
8. **reference_aa** — Reference amino acid
9. **variant_aa** — Variant amino acid
10. **notes** — Additional annotations

## Validation

Verify the database is properly indexed:

```bash
# Check tabix index
tabix -H data/herpesdrg_cmv.bed.gz

# Query specific region (UL97 kinase domain)
tabix data/herpesdrg_cmv.bed.gz NC_006273:94950-95200
```

## References

1. Tilloy V, et al. (2024) "Comprehensive Herpesviruses Antiviral drug Resistance Mutation Database (CHARMD)." *Antiviral Research* 231:106016. PMID: 39349222
2. Charles O, et al. HerpesDRG Database. https://github.com/ojcharles/herpesdrg-db
3. Chou S, et al. (2021) "Cytomegalovirus drug resistance: a review." *Antiviral Therapy*
4. CNR Herpesvirus, CHU Limoges. CODEX MV / CHARMD. https://www.unilim.fr/cnr-herpesvirus/outils/codexmv/

## License

Data derived from HerpesDRG (MIT License per GitHub repository).

## Contact

For questions about this database, open an issue on GitHub.
