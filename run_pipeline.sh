#!/bin/bash
# CMV HerpesDRG vcfanno pipeline
# Annotate CMV whole-genome VCFs with drug resistance mutations
# Uses NC_006273.2 (current AD169 reference) coordinates
set -e

# Build database from HerpesDRG TSV
python3 scripts/build_database.py \
    --herpesdrg data/herpesdrg-db.tsv \
    --genome data/AD169_NC_006273.2.fa \
    --genbank data/AD169_NC_006273.2.gb \
    --output data/herpesdrg_cmv.bed

# Create ref/alt version (adds REF/ALT columns for vcfanno base-level matching)
python3 scripts/add_ref_alt_to_bed.py \
    --input data/herpesdrg_cmv.bed \
    --genome data/AD169_NC_006273.2.fa \
    --genbank data/AD169_NC_006273.2.gb \
    --output data/herpesdrg_cmv_refalt.bed

# Compress and tabix-index the ref/alt BED (single bgz block + tbi)
python3 - <<'PY'
import pysam, os, gzip, struct
base = "data/herpesdrg_cmv_refalt.bed"
gz = base + ".gz"
tbi = gz + ".tbi"
for f in [gz, tbi]:
    if os.path.exists(f): os.remove(f)
pysam.tabix_compress(base, gz, force=True)
pysam.tabix_index(gz, preset="bed", force=True)
with gzip.open(gz, "rt") as f:
    starts = [int(l.split("\t")[1]) for l in f if not l.startswith("#")]
print(f"Bed: sorted={starts==sorted(starts)}, unique={len(set(starts))}/{len(starts)}")
PY

# Annotate VCF (uses ref/alt BED for base-level matching)
vcfanno -p 1 data/vcfanno_cmv_refalt.toml input.vcf.gz > output.vcf
