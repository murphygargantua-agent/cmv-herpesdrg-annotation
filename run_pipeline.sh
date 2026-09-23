#!/bin/bash
# CMV HerpesDRG vcfanno annotation pipeline
# Usage: ./run_pipeline.sh <herpesdrg.tsv> <ad169.fa> <ad169.gb> <input.vcf.gz> <output.vcf.gz>
set -e

HERPESDRG=${1:-data/herpesdrg-db.tsv}
GENOME=${2:-/tmp/AD169.fa}
GENBANK=${3:-/tmp/AD169_X17403.gb}
INPUT_VCF=${4:-input.vcf.gz}
OUTPUT_VCF=${5:-annotated.vcf.gz}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"

echo "=== Step 1: Build BED ==="
python3 "${SCRIPT_DIR}/scripts/build_database.py" \
    --herpesdrg "$HERPESDRG" \
    --genome "$GENOME" \
    --genbank "$GENBANK" \
    --output "${DATA_DIR}/herpesdrg_cmv.bed"

echo "=== Step 2: Deduplicate, sort, bgzip, index ==="
python3 - <<PY
import pysam, os
from collections import defaultdict
bed = "${DATA_DIR}/herpesdrg_cmv.bed"
merged = defaultdict(list)
with open(bed) as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.strip("\n").split("\t")
        merged[(p[0], int(p[1]), int(p[2]))].append(p[3])
records = sorted(merged.items(), key=lambda x: (x[0][0], x[0][1]))
out = bed + ".merged"
with open(out, "w") as f:
    f.write("##bedFormat=4\n")
    for (c, s, e), labels in records:
        f.write(f"{c}\t{s}\t{e}\t{';'.join(labels)}\n")
bgz = "${DATA_DIR}/herpesdrg_cmv.bed.gz"
for p in [bgz, bgz + ".tbi"]:
    if os.path.exists(p): os.remove(p)
pysam.tabix_compress(out, bgz, force=True)
pysam.tabix_index(bgz, seq_col=0, start_col=1, end_col=2, force=True)
print(f"Wrote {len(records)} unique positions to {bgz}")
PY

echo "=== Step 3: Annotate VCF with vcfanno ==="
TOML="${DATA_DIR}/vcfanno_cmv.toml"
/workspace/tools/vcfanno -p 1 "$TOML" "$INPUT_VCF" > "${OUTPUT_VCF%.gz}"

echo "=== Done ==="
echo "Output: ${OUTPUT_VCF}"
