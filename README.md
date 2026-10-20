# CMV HerpesDRG vcfanno Database

A database of CMV antiviral resistance mutations from HerpesDRG, mapped to AD169/NC_006273.2 coordinates, formatted for vcfanno.

## Quick Start

```bash
# Build the BED (requires AD169 reference files)
python3 scripts/build_database.py \
    --herpesdrg data/herpesdrg-db.tsv \
    --genome data/AD169_NC_006273.2.fa \
    --genbank data/AD169_NC_006273.2.gb \
    --output data/herpesdrg_cmv.bed

# Create ref/alt version
python3 scripts/add_ref_alt_to_bed.py \
    --input data/herpesdrg_cmv.bed \
    --genome data/AD169_NC_006273.2.fa \
    --genbank data/AD169_NC_006273.2.gb \
    --output data/herpesdrg_cmv_refalt.bed

# Sort + index
bash run_pipeline.sh
# or: INPUT=your_variants.vcf.gz bash run_pipeline.sh
```

## Database Contents

| Metric | Value |
|---|---|
| Total mutations | 772 (mapped to 386 unique codon positions) |
| Genes | UL54, UL97, UL56, UL51, UL27, UL89 |
| Reference | AD169 (NC_006273.2, 235,646 bp) |
| Coverage | ~100% (386 unique positions) |

## vcfanno Config

```toml
[[annotation]]
file = "data/herpesdrg_cmv_refalt.bed.gz"
columns = [4]
names = ["cmv_resistance"]
ops = ["uniq"]
type = "String"
```

## Output Format

`UL54_D301N|Ganciclovir=2.6;Cidofovir=3;Foscarnet=0.5|UL54_D301N|Aciclovir=0.2;Cidofovir=14`

## Critical Requirements

- **BED must be sorted and tabix-indexed** (bgzip + tabix, preset=bed) — `run_pipeline.sh` sorts automatically
- **VCF must be sorted by position** (numerically, 5'→3') — vcfanno panics if out of order; `run_pipeline.sh` sorts automatically via `bcftools sort`
- **Chromosome names must match** — this DB uses `NC_006273.2`
- **Deduplicate intervals** — `build_database.py` merges duplicates

## Known Issues

- **Amino acid mismatch warnings**: ~4.4% of mutations (34/772) have reference amino acids that differ between the NC_006273.2 (Merlin) reference and the HerpesDRG AD169 reference. These are flagged with a warning but still annotated. The mismatches are concentrated in UL89.
- **HerpesDRG BED has duplicate genomic positions**: `ops = ["uniq"]` keeps only the first match.
- **minMutFinder vs vcfanno**: minMutFinder maps to CDS reference and calls variants de novo; vcfanno annotates pre-called VCF variants against the BED database. They use different approaches but the same reference coordinates.

## Installation

vcfanno binary: https://github.com/brentp/vcfanno/releases/download/v0.3.9/vcfanno_0.3.9_linux_amd64.zip

## How It Works

1. `build_database.py` parses HerpesDRG TSV (772 target gene rows) and AD169 GenBank annotation, computes codon-spanning genomic intervals, deduplicates entries, writes a sorted BED file.
2. `add_ref_alt_to_bed.py` reads the BED, looks up REF/ALT bases from the AD169 reference, writes a BED with REF/ALT columns for base-level vcfanno matching.
3. `run_pipeline.sh` sorts and tabix-indexes the BED, sorts the VCF, and runs vcfanno.
4. vcfanno annotates VCF variants that overlap the BED intervals.
