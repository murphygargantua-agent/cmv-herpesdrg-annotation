#!/bin/bash
# CMV HerpesDRG vcfanno pipeline
# Annotate CMV whole-genome VCFs with drug resistance mutations
# Uses NC_006273.2 (current AD169 reference) coordinates
set -e

# Ensure vcfanno is available (installed in /root/miniconda3/bin)
# and minmutfinder env Python (with pysam) takes precedence
export PATH="/root/miniconda3/envs/minmutfinder/bin:/root/miniconda3/bin:$PATH"

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

# Sort BED files (safety check) and compress + tabix-index
python3 - <<'PY'
import pysam, os, gzip
for base in ["data/herpesdrg_cmv.bed", "data/herpesdrg_cmv_refalt.bed"]:
    with open(base) as f:
        lines = f.readlines()
    header = [l for l in lines if l.startswith("#")]
    data = [l for l in lines if not l.startswith("#")]
    data.sort(key=lambda x: (x.split("\t")[0], int(x.split("\t")[1])))
    with open(base + ".sorted", "w") as f:
        f.writelines(header + data)
    gz = base + ".sorted.gz"
    tbi = gz + ".tbi"
    for f in [gz, tbi]:
        if os.path.exists(f): os.remove(f)
    pysam.tabix_compress(base + ".sorted", gz, force=True)
    pysam.tabix_index(gz, preset="bed", force=True)
    with gzip.open(gz, "rt") as f:
        starts = [int(l.split("\t")[1]) for l in f if not l.startswith("#")]
    print(f"{base}: sorted={starts==sorted(starts)}, unique={len(set(starts))}/{len(starts)}")
PY

# Ensure input VCF is bgzipped
if [[ "$INPUT" == *.vcf.gz ]]; then
    INPUT_GZ="$INPUT"
elif [[ "$INPUT" == *.vcf ]]; then
    INPUT_GZ="${INPUT}.gz"
    bcftools view -Oz -o "$INPUT_GZ" "$INPUT"
else
    echo "Error: INPUT must be .vcf or .vcf.gz" >&2
    exit 1
fi

# Sort VCF by position (required by vcfanno — will panic if unsorted)
SORTED_VCF="${INPUT%.vcf.gz}.sorted.vcf.gz"
bcftools sort -Oz -o "$SORTED_VCF" "$INPUT_GZ"
tabix -p vcf "$SORTED_VCF"

# Annotate VCF with ref/alt BED (base-level matching)
# -permissive-overlap: annotate even when REF/ALT differs from BED (handles 7.8% AD169 vs Merlin strain differences)
vcfanno -p 1 -permissive-overlap data/vcfanno_cmv_refalt.toml "$SORTED_VCF" > output.vcf
