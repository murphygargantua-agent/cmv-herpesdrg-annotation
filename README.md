# CMV HerpesDRG vcfanno Database

A database of CMV antiviral resistance mutations from HerpesDRG, mapped to AD169/NC_006273.2 coordinates, formatted for `vcfanno`.

## Quick Start

```bash
# Build the BED (requires AD169 reference files)
python3 scripts/build_database.py \
  --herpesdrg herpesdrg-db.tsv \
  --genome AD169.fasta \
  --genbank AD169_NC_006273.gb \
  --output data/herpesdrg_cmv.bed

# Index for vcfanno (bgzip + tabix)
python3 scripts/index_bed.py

# Annotate
vcfanno -p 1 data/vcfanno_cmv.toml your_variants.vcf.gz > annotated.vcf
```

## Database Contents

| Metric | Value |
|--------|-------|
| Total mutations | 749 (mapped to 373 unique codon positions) |
| Genes | UL54, UL97, UL56, UL51, UL27, UL89 |
| Reference | AD169 (NC_006273.2, 235,646 bp) |
| Coverage | ~100% (373 unique positions) |

## vcfanno Config (data/vcfanno_cmv.toml)

```toml
[[annotation]]
file = "data/herpesdrg_cmv.bed.gz"
columns = [4]
names = ["cmv_resistance"]
ops = ["uniq"]
type = "String"
```

## Output Format

```
UL54_D301N|Ganciclovir=2.6;Cidofovir=3;Foscarnet=0.5|UL54_D301N|Aciclovir=0.2;Cidofovir=14
```

## Critical Requirements

1. **BED must be sorted and tabix-indexed** (bgzip + tabix, `preset=bed`)
2. **Deduplicate intervals** at the same start position — `build_database.py` merges duplicates
3. **Use `columns = [4]` and `ops = ["uniq"]`** in the TOML config (not `fields`)
4. **VCF must be sorted by position** (numerically, 5'→3') — vcfanno fails if out of order
5. **Chromosome names must match** — this DB uses `NC_006273.2`

## Installation

vcfanno binary: `https://github.com/brentp/vcfanno/releases/download/v0.3.9/vcfanno_0.3.9_linux_amd64.zip`

## How It Works

1. `build_database.py` parses HerpesDRG TSV (779 target gene rows) and AD169 GenBank annotation
2. Computes codon-spanning 3-base genomic intervals using the CDS coordinates
3. Deduplicates entries at the same genomic position (vcfanno fails with duplicate starts)
4. Writes a BED file with `bedFormat=4` comment
5. `index_bed.py` bgzips and tabix-indexes the BED file
6. vcfanno annotates VCF variants that overlap the BED intervals

## Comparison with CHARMD

| Feature | HerpesDRG | CHARMD |
|---------|-----------|--------|
| Source | GitHub TSV (ojcharles/herpesdrg-db) | CHARMD (CNR Herpesvirus, France) |
| Mutations | 749 (373 unique) | 441 (355 unique) |
| Coverage | ~100% | 71.3% |
| Drugs | 8 antivirals | GCV, MBV, FOS, CDD, BCD |
| EC50 data | Yes | Yes |
| Phenotype | Yes | Yes |

Both databases use the same AD169 reference (NC_006273.2) and can be combined by concatenating the BED files (after deduplication).
