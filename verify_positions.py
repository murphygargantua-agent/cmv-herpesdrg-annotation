#!/usr/bin/env python3
"""
verify_positions.py - Verify that vcfanno BED positions are correct for NC_006273.2.

Checks that the amino acid positions in the HerpesDRG BED map to the correct
genomic coordinates in the NC_006273.2 (Merlin/AD169) reference, using the
CDS coordinates from the GenBank file.

Usage:
    python verify_positions.py --bed data/herpesdrg_cmv_refalt.bed \
        --genome data/AD169_NC_006273.2.fa --genbank data/AD169_NC_006273.2.gb
"""
import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", required=True, help="Path to BED file")
    parser.add_argument("--genome", required=True, help="Path to FASTA reference")
    parser.add_argument("--genbank", required=True, help="Path to GenBank file")
    args = parser.parse_args()

    ref = next(SeqIO.parse(args.genome, "fasta")).seq.upper()
    gb = SeqIO.read(args.genbank, "genbank")

    cds = {}
    for feat in gb.features:
        if feat.type == "CDS":
            gene = feat.qualifiers.get("gene", [""])[0]
            cds[gene] = (int(feat.location.start), int(feat.location.end), feat.location.strand)

    mismatches = []
    total = 0
    with open(args.bed) as f:
        f.readline()  # skip header line
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 6:
                continue
            chrom, start, end, ann, ref_base, alt_base = parts[:6]
            mut_name = ann.split("|")[0]
            if "_" not in mut_name:
                continue
            gene = mut_name.split("_")[0]
            aa_pos = int("".join(filter(str.isdigit, mut_name.split("_")[1])))
            total += 1
            if gene not in cds:
                continue
            cds_start, cds_end, strand = cds[gene]
            if strand == 1:
                codon_start = cds_start + (aa_pos - 1) * 3
            else:
                codon_start = cds_end - aa_pos * 3
            if codon_start + 3 > len(ref):
                continue
            ref_bases = ref[codon_start:codon_start + 3]
            ref_codon = ref_bases if strand == 1 else str(Seq(ref_bases).reverse_complement())
            ref_aa = str(Seq(ref_codon).translate())
            if ref_aa != mut_name.split("_")[1][0]:
                mismatches.append((mut_name, aa_pos, ref_codon, ref_aa, ref_base))

    print(f"Total: {total}, Mismatches: {len(mismatches)} ({len(mismatches)/total*100:.1f}%)")
    if mismatches:
        print("\nReference AA mismatches (HerpesDRG vs NC_006273.2):")
        for m in mismatches[:10]:
            print(f"  {m[0]}: codon={m[2]}, translated={m[3]}, BED_REF={m[4]}")


if __name__ == "__main__":
    main()