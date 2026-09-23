# CMV Antiviral Resistance Variant Database

A curated database of cytomegalovirus (CMV) antiviral resistance mutations
formatted for integration with `vcfanno`. Built from [HerpesDRG](https://github.com/ojcharles/herpesdrg-db)
data, mapped to the [AD169](https://www.ncbi.nlm.nih.gov/nuccore/X17403) reference genome
(NCBI accession X17403.1, NC_006273).

## Quick Start

### Install vcfanno

```bash
# Download the static binary from GitHub releases:
curl -L -o /workspace/tools/vcfanno https://github.com/brentp/vcfanno/releases/download/v0.3.9/vcfanno_0.3.9_linux_amd64.zip
unzip vcfanno_0.3.9_linux_amd64.zip
chmod +x vcfanno
```

### Build the database

```bash
python3 scripts/build_database.py \
    --herpesdrg herpesdrg-db.tsv \
    --genome AD169.fasta \
    --genbank AD169_annotation.gb \
    --output data/herpesdrg_cmv.bed
```

The script:
1. Parses AD169 GenBank CDS annotations for target genes (UL54, UL97, UL56, UL51, UL27, UL89).
2. Derives amino-acid → genomic codon coordinates using the `aa_to_genomic` mapping.
3. Verifies each codon translates to the expected reference amino acid against the AD169 reference.
4. Maps each mutation to its 3-base genomic codon interval and attaches drug/EC50 data.
5. Deduplicates entries at the same genomic position.

### Compress and index

```bash
# Merge duplicates, sort, bgzip, and tabix-index:
python3 - <<'PY'
import pysam, os
from collections import defaultdict
merged = defaultdict(list)
with open("data/herpesdrg_cmv.bed") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.strip().split("\t")
        merged[(p[0], int(p[1]), int(p[2]))].append(p[3])
records = sorted(merged.items(), key=lambda x: (x[0][0], x[0][1]))
with open("data/herpesdrg_cmv_merged.bed", "w") as f:
    f.write("##bedFormat=4\n")
    for (c, s, e), labels in records:
        f.write(f"{c}\t{s}\t{e}\t{';'.join(labels)}\n")
pysam.tabix_compress("data/herpesdrg_cmv_merged.bed", "data/herpesdrg_cmv.bed.gz", force=True)
pysam.tabix_index("data/herpesdrg_cmv.bed.gz", seq_col=0, start_col=1, end_col=2, force=True)
PY
```

### Annotate a VCF

```bash
vcfanno -p 1 data/vcfanno_cmv.toml your_variants.vcf.gz > annotated.vcf
```

**vcfanno_cmv.toml:**
```toml
[[annotation]]
file = "data/herpesdrg_cmv.bed.gz"
columns = [4]
names = ["cmv_resistance"]
ops = ["uniq"]
```

**Important:** The VCF must be sorted by position (numerically, 5'→3'). vcfanno
will fail with "intervals out of order" if the VCF is out of order.

### Output Format

Annotations appear in the VCF INFO field as:
```
cmv_resistance=UL54_D301N|Ganciclovir=2.6,Cidofovir=3,Foscarnet=0.5
```

### Expected VCF format

```vcf
##fileformat=VCFv4.2
##contig=<ID=NC_006273,length=229354>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
NC_006273	902	.	G	A	.	.	.	.	.
NC_006273	79729	.	T	C	.	.	.	.	.
```

## Database Contents

| Metric | Value |
|--------|-------|
| Total mutations in DB | 716 (verified against AD169) |
| Unique codon positions | 345 (after dedup) |
| Genes covered | UL97, UL54, UL27, UL56, UL51, UL89 |
| Reference genome | AD169 (NC_006273 / X17403.1) |
| Source | HerpesDRG (LANE-CMR) |
| Coverage | 97.3% (716/735 mapped rows) |

### Gene Distribution

| Gene | Count | Drug Classes |
|------|-------|--------------|
| UL54 | ~270 | Ganciclovir, Foscarnet, Cidofovir, Brincidofovir, Valganciclovir |
| UL97 | ~184 | Ganciclovir, Valganciclovir, Maribavir |
| UL56 | ~104 | Letermovir |
| UL51 | ~9 | Letermovir |
| UL27 | ~19 | Maribavir |
| UL89 | — | Letermovir |

## Data Sources

### Primary: HerpesDRG (LANE-CMR)

- **Source:** [HerpesDRG Database](https://github.com/ojcharles/herpesdrg-db)
- **Download:** `herpesdrg-db.tsv`
- **Coverage:** 602 HCMV mutations (586 mapped = 97.3%)

### Reference: AD169 (NC_006273 / X17403.1)

- **Genome FASTA:** `NCBI Entrez EFetch db=nucleotide id=X17403.1 rettype=fasta`
- **GenBank annotation:** `NCBI Entrez EFetch db=nucleotide id=X17403.1 rettype=gb`

## Known Limitations

1. **v0.3.9 issue:** vcfanno 0.3.9 may fail with "intervals out of order" if the BED file
   has duplicate start positions — always deduplicate/merge first.
2. **REF/ALT matching:** By default vcfanno requires exact REF/ALT match between VCF and BED.
   Use `-permissive-overlap` to annotate without matching alleles.
3. **Strand handling:** For reverse-strand genes (UL54, UL97, UL56, UL51, UL27, UL89),
   the codon is reverse-complemented; verification catches any coordinate errors.
4. **Strain differences:** AD169-specific coordinates; other strains (Merlin, Toledo)
   require separate coordinate mapping.
5. **Unmapped mutations:** 33 mutations failed codon verification (likely due to
   reference amino acid mismatches or positions beyond AD169 protein length).
6. **No VCF→DRM classification:** vcfanno only labels which known resistance mutation
   a variant corresponds to; it does not classify phenotype severity. Use CHARMD or
   genotype-phenotype tables for that.
