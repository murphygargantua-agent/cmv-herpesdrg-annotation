#!/usr/bin/env python3
"""
build_database.py - Rebuild the CMV HerpesDRG vcfanno database.

Usage:
    python build_database.py --herpesdrg herpesdrg-db.tsv \
                             --genome AD169.fasta \
                             --genbank AD169_annotation.gb \
                             --output herpesdrg_cmv.bed
"""
import argparse
import csv
import json
import sys
from pathlib import Path


def parse_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                seq.append(line)
    return ''.join(seq)


def parse_genbank_genes(path, target_genes):
    genes = {}
    try:
        from Bio import GenBank
        with open(path) as f:
            record = GenBank.read(f)
        for feature in record.features:
            if feature.type == "CDS" and "gene" in feature.qualifiers:
                raw = feature.qualifiers["gene"][0]
                if raw.startswith("HCMV"):
                    raw = raw[4:]
                if raw in target_genes:
                    loc = feature.location
                    genes[raw] = {
                        "start": loc.start + 1,
                        "end": loc.end,
                        "length": loc.end - loc.start,
                    }
    except ImportError:
        json_path = Path(path).parent / "ad169_gene_mapping.json"
        if json_path.exists():
            with open(json_path) as f:
                genes = json.load(f)
    return genes


def parse_mutation(aa_change):
    if 'del' in aa_change.lower():
        part = aa_change.replace("del", "").strip()
        start, end = part.split("-")
        return f"del_{start}_{end}", int(start), None
    ref = aa_change[0]
    end = aa_change[-1]
    pos_str = aa_change[1:-1]
    pos = int(pos_str)
    return ref, pos, end


def aa_to_nt(gene, aa_pos, gene_coords, genome_seq):
    coords = gene_coords[gene]
    gene_start = coords["start"]
    nt_start = (gene_start - 1) + (aa_pos * 3)
    ref_nt = genome_seq[nt_start:nt_start + 3]
    return nt_start, nt_start + 3, ref_nt


def build_database(herpesdrg_path, genome_path, genbank_path, output_path, target_genes):
    print(f"Loading genome from {genome_path}...")
    genome_seq = parse_fasta(genome_path)
    print(f"  Length: {len(genome_seq)} bp")

    print(f"Loading gene coordinates from {genbank_path}...")
    gene_coords = parse_genbank_genes(genbank_path, target_genes)
    print(f"  Found {len(gene_coords)} genes: {', '.join(gene_coords)}")

    print(f"Loading HerpesDRG from {herpesdrg_path}...")
    mutations = []
    with open(herpesdrg_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get("virus", "").strip().upper() != "HCMV":
                continue
            if row.get("status", "").strip().upper() != "A":
                continue
            if row.get("gene", "").strip() not in target_genes:
                continue
            mutations.append(row)
    print(f"  {len(mutations)} active HCMV mutations in target genes")

    print("Mapping mutations to genomic coordinates...")
    mapped = []
    unmapped = []

    for mut in mutations:
        gene = mut["gene"].strip()
        aa_change = mut["amino_acid_change"].strip()
        ref_aa, pos, var_aa = parse_mutation(aa_change)

        if pos is None:
            unmapped.append(f"{gene}_{aa_change} (deletion/unparseable)")
            continue

        try:
            nt_start, nt_end, ref_nt = aa_to_nt(gene, pos, gene_coords, genome_seq)

            bed_line = [
                "NC_006273",
                nt_start,
                nt_end,
                f"{gene}_{aa_change}",
                gene,
                mut.get("drug", "").strip(),
                mut.get("fold_change", "N/A").strip(),
                ref_aa,
                var_aa or "-",
                mut.get("notes", "").strip(),
            ]
            mapped.append(bed_line)

        except (KeyError, IndexError, ValueError) as e:
            unmapped.append(f"{gene}_{aa_change} ({e})")

    print(f"Writing BED file to {output_path}...")
    print(f"  Mapped: {len(mapped)}")
    print(f"  Unmapped: {len(unmapped)}")

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter='\t')
        for line in mapped:
            writer.writerow(line)

    if unmapped:
        unmapped_path = Path(output_path).with_suffix(".bed.unmapped")
        with open(unmapped_path, "w") as f:
            for m in unmapped:
                f.write(m + "\n")
        print(f"  Unmapped list saved to {unmapped_path}")

    return len(mapped), len(unmapped)


def main():
    parser = argparse.ArgumentParser(
        description="Build CMV HerpesDRG vcfanno database"
    )
    parser.add_argument(
        "--herpesdrg", required=True,
        help="Path to herpesdrg-db.tsv"
    )
    parser.add_argument(
        "--genome", required=True,
        help="Path to AD169 FASTA"
    )
    parser.add_argument(
        "--genbank", required=True,
        help="Path to AD169 GenBank"
    )
    parser.add_argument(
        "--output", default="herpesdrg_cmv.bed",
        help="Output BED file"
    )
    args = parser.parse_args()

    target_genes = ["UL97", "UL54", "UL27", "UL56", "UL51", "UL89"]
    mapped, unmapped = build_database(
args.herp
        args.genome,
        args.genbank,
        args.output,
        target_genes,
    )
    print(f"Done! Mapped: {mapped}, Unmapped: {unmapped}")
