#!/usr/bin/env python3
"""
add_ref_alt.py - Add REF/ALT columns to vcfanno BED files.

vcfanno supports matching on REF/ALT bases when the BED header contains
'##chr ref alt' columns. This script:
  1. Reads the existing BED file (with codon-spanning intervals)
  2. Looks up the 3 reference bases from the NC_006273.2 AD169 reference
  3. Computes the mutation bases (alt) from the codon label
  4. Writes a new BED file with REF/ALT columns

Usage:
    python add_ref_alt.py --input herpesdrg_cmv.bed --genome AD169_NC_006273.2.fa \
                         --genbank AD169_NC_006273.2.gb \
                         --output herpesdrg_cmv_refalt.bed
"""

import argparse
import re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


def load_genome(fasta_path):
    record = SeqIO.read(fasta_path, "fasta")
    return str(record.seq).upper()


def codon_to_aa(codon):
    return str(Seq(codon.upper()).translate())


def parse_label(label):
    """Parse 'UL54_D301N|drug=val' -> (gene, ref_aa, pos, mut_aa, extra)."""
    main = label.split("|")[0]
    m = re.match(r"^([A-Z]+[0-9]+)_(\D+)(\d+)(\D+)$", main)
    if not m:
        return None
    return (m.group(1), m.group(2), int(m.group(3)), m.group(4),
            "|".join(label.split("|")[1:]) if "|" in label else "")


def build_codon_table(genome, genbank_path, target_genes):
    """Build codon table from GenBank: gene -> {strand, codons{pos: genomic_start}}."""
    gb = SeqIO.read(genbank_path, "genbank")
    codon_table = {}
    for feat in gb.features:
        if feat.type != "CDS":
            continue
        gene_q = feat.qualifiers.get("gene", [""])[0]
        gene_name = gene_q[4:].split()[0] if gene_q.startswith("HCMV") else gene_q.split()[0]
        if gene_name not in target_genes:
            continue
        cds_start = int(feat.location.start)
        cds_end = int(feat.location.end)
        strand = feat.location.strand
        num_codons = (cds_end - cds_start) // 3
        table = {}
        if strand == 1:
            for n in range(1, num_codons + 1):
                table[n] = cds_start + (n - 1) * 3
        else:
            for n in range(1, num_codons + 1):
                table[n] = cds_end - n * 3
        codon_table[gene_name] = {"strand": strand, "codons": table}
    return codon_table


def find_alt_codon(ref_codon, mut_aa):
    """Find codon that translates to mut_aa, differs from ref_codon by 1 nt."""
    codon_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    candidates = [c for c in codon_table if codon_table[c] == mut_aa]
    for c in candidates:
        if c == ref_codon:
            continue
        if sum(1 for a, b in zip(ref_codon, c) if a != b) == 1:
            return c
    return candidates[0] if candidates else ref_codon


def main():
    parser = argparse.ArgumentParser(
        description="Add REF/ALT columns to vcfanno BED database")
    parser.add_argument("--input", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--genbank", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--target-genes", nargs="+",
                        default=["UL54", "UL97", "UL56", "UL51", "UL27", "UL89"])
    args = parser.parse_args()

    genome = load_genome(args.genome)
    codon_table = build_codon_table(genome, args.genbank, args.target_genes)
    print(f"Genes: {sorted(codon_table.keys())}")

    mutations = []
    headers = []
    with open(args.input) as f:
        for line in f:
            if line.startswith("##"):
                headers.append(line)
                continue
            parts = line.strip().split("\t")
            mutations.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
    print(f"Total entries: {len(mutations)}")

    added = failed = 0
    out = []
    for chrom, start, end, label in mutations:
        gene_m = re.match(r"^([A-Z]+[0-9]+)", label.split("|")[0])
        gene = gene_m.group(1) if gene_m else None

        if not gene or gene not in codon_table:
            out.append(f"{chrom}\t{start}\t{end}\t{label}\t.\t.")
            failed += 1
            continue

        info = codon_table[gene]
        strand = info["strand"]
        codons = info["codons"]

        codon_num = None
        for n, pos in codons.items():
            if pos == start:
                codon_num = n
                break
        if codon_num is None:
            out.append(f"{chrom}\t{start}\t{end}\t{label}\t.\t.")
            failed += 1
            continue

        if start + 3 > len(genome):
            out.append(f"{chrom}\t{start}\t{end}\t{label}\t.\t.")
            failed += 1
            continue

        ref_bases = genome[start:start + 3]
        parsed = parse_label(label)
        if parsed is None:
            out.append(f"{chrom}\t{start}\t{end}\t{label}\t.\t.")
            failed += 1
            continue

        gene_p, ref_aa, pos, mut_aa, extra = parsed
        if pos != codon_num:
            out.append(f"{chrom}\t{start}\t{end}\t{label}\t.\t.")
            failed += 1
            continue

        ref_codon = ref_bases if strand == 1 else str(Seq(ref_bases).reverse_complement())
        ref_aa_check = codon_to_aa(ref_codon)
        if ref_aa_check != ref_aa:
            print(f"WARNING: {label}: ref_aa {ref_aa} != translated {ref_aa_check}")

        alt_codon = find_alt_codon(ref_codon, mut_aa)
        alt_genomic = alt_codon.upper() if strand == 1 else str(Seq(alt_codon).reverse_complement()).upper()

        out.append(f"{chrom}\t{start}\t{end}\t{label}\t{ref_bases.upper()}\t{alt_genomic}")
        added += 1

    with open(args.output, "w") as f:
        f.write("##chr\tref\talt\n")
        for h in headers:
            f.write(h)
        for line in out:
            f.write(line + "\n")

    print(f"Added: {added}, Failed: {failed}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()