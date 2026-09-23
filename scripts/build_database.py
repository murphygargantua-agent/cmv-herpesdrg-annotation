#!/usr/bin/env python3
"""
build_database.py - Rebuild the CMV HerpesDRG vcfanno database.

Usage:
    python build_database.py --herpesdrg herpesdrg-db.tsv \\
                             --genome AD169.fasta \\
                             --genbank AD169_annotation.gb \\
                             --output herpesdrg_cmv.bed

Downloads:
    - HerpesDRG: https://github.com/ojcharles/herpesdrg-db/raw/main/herpesdrg-db.tsv
    - AD169 genome: NCBI accession X17403.1 (via Entrez EFetch)
    - AD169 annotation: NCBI accession AE016761 (via Entrez EFetch)
"""

import argparse
import csv
import re
import os
import sys
from pathlib import Path

from Bio import Entrez, SeqIO


def parse_fasta(path: str) -> str:
    """Extract genome sequence from FASTA file."""
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                seq.append(line)
    return ''.join(seq)


def parse_ad169_coordinates(target_genes: list[str], genbank_path: str) -> dict:
    """
    Parse AD169 GenBank annotation to get CDS codon coordinates.
    Returns dict mapping gene name to codon-start positions.
    For forward strand: codon n starts at cds_start + (n-1)*3
    For reverse strand: codon n starts at cds_end - (n-1)*3
    """
    from Bio.Seq import Seq
    with open(genbank_path) as f:
        record = SeqIO.read(f, 'genbank')

    gene_codon_coords = {}
    for feature in record.features:
        if feature.type == "CDS" and "gene" in feature.qualifiers:
            raw_gene = feature.qualifiers["gene"][0]
            if raw_gene.startswith("HCMV"):
                gene_name = raw_gene[4:].split()[0]
            else:
                gene_name = raw_gene.split()[0]

            if gene_name not in target_genes:
                continue

            cds_start = int(feature.location.start)  # 0-based
            cds_end = int(feature.location.end)      # 1-based
            strand = feature.location.strand
            num_codons = (cds_end - cds_start) // 3

            # Build dict: aa_pos -> codon_genomic_start (0-based, first base of codon)
            codon_map = {}
            if strand == 1:
                # Forward strand: codon n starts at cds_start + (n-1)*3
                for n in range(1, num_codons + 1):
                    codon_map[n] = cds_start + (n - 1) * 3
            else:
                # Reverse strand: codon n starts at cds_end - n*3 (then we use the start as 0-based)
                for n in range(1, num_codons + 1):
                    # Codon n (in protein order) on reverse strand occupies [cds_end - n*3, cds_end - (n-1)*3)
                    codon_map[n] = cds_end - n * 3

            gene_codon_coords[gene_name] = codon_map

    return gene_codon_coords


def build_database(herpesdrg_path: str, genome_path: str,
                   genbank_path: str, output_path: str,
                   target_genes: list[str]) -> tuple[int, int]:
    """
    Build CMV mutation database from HerpesDRG TSV.

    Args:
        herpesdrg_path: Path to herpesdrg-db.tsv
        genome_path: Path to AD169 FASTA
        genbank_path: Path to AD169 GenBank annotation
        output_path: Output BED file path
        target_genes: List of gene names to include

    Returns:
        (mapped_count, unmapped_count)
    """
    # Parse genome
    genome_seq = parse_fasta(genome_path)
    print(f"Genome: {len(genome_seq)} bp")

    # Parse gene coordinates
    print("Parsing GenBank annotations...")
    gene_codon_coords = parse_ad169_coordinates(target_genes, genbank_path)
    print(f"Found {len(gene_codon_coords)} genes:")
    for gene in sorted(gene_codon_coords.keys()):
        codons = gene_codon_coords[gene]
        n_codons = len(codons)
        print(f"  {gene}: {n_codons} codons")

    # Parse HerpesDRG TSV
    print("\nLoading HerpesDRG TSV...")
    with open(herpesdrg_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)

    print(f"Total rows: {len(rows)}")
    target_rows = [r for r in rows if r['gene'] in target_genes]
    print(f"Target gene rows: {len(target_rows)}")

    # Build mutation list
    print("\nBuilding mutation list...")
    cmv_mutations = []
    unmapped = 0
    unparsed = 0

    for row in target_rows:
        gene = row['gene']
        aa_change = row['aa_change']

        # Parse aa_change format: RefPosMut (e.g., D301N = D at pos 301 -> N)
        match = re.match(r'^([A-Z])(\d+)([A-Z])$', aa_change)
        if not match:
            unparsed += 1
            continue

        ref_aa = match.group(1)
        pos = int(match.group(2))
        mut_aa = match.group(3)

        if gene not in gene_codon_coords:
            unmapped += 1
            continue

        coords = gene_codon_coords[gene]
        if pos not in coords:
            unmapped += 1
            continue

        genomic_pos = coords[pos]  # 0-based first base of codon

        # Build drug resistance info string
        drug_info = []
        for drug in ['Ganciclovir', 'Aciclovir', 'Cidofovir', 'Foscarnet',
                     'Brincidofovir', 'Letermovir', 'Maribavir']:
            val = row.get(drug, '')
            if val and str(val).strip():
                drug_info.append(f"{drug}={val}")

        # Label format: GENE_AA_CHANGE|Drug1=val1;Drug2=val2
        label = f"{gene}_{aa_change}"
        if drug_info:
            label += '|' + ';'.join(drug_info)

        chrom = 'NC_006273'  # AD169 NCBI accession
        cmv_mutations.append({
            'chrom': chrom,
            'start': genomic_pos,
            'end': genomic_pos + 1,
            'label': label
        })

    # Write BED file, sorted and merged by genomic position
    from collections import defaultdict
    merged = defaultdict(list)
    for m in cmv_mutations:
        key = (m['chrom'], m['start'], m['end'])
        merged[key].append(m['label'])

    # Sort by chrom, start, end
    sorted_keys = sorted(merged.keys(), key=lambda x: (x[0], x[1], x[2]))

    print(f"\nMapped: {len(cmv_mutations)}")
    print(f"Unparsed aa_change: {unparsed}")
    print(f"Unmapped: {unmapped}")
    print(f"Unique genomic positions: {len(merged)}")

    with open(output_path, 'w') as f:
        f.write("##bedFormat=4\n")
        for key in sorted_keys:
            chrom, start, end = key
            # Deduplicate labels
            unique_labels = []
            seen = set()
            for label in merged[key]:
                if label not in seen:
                    seen.add(label)
                    unique_labels.append(label)
            f.write(f"{chrom}\t{start}\t{end}\t{'|'.join(unique_labels)}\n")

    print(f"\nWrote {len(sorted_keys)} unique positions to {output_path}")

    return len(cmv_mutations), unmapped + unparsed


def main():
    parser = argparse.ArgumentParser(
        description="Build CMV HerpesDRG vcfanno database"
    )
    parser.add_argument(
        '--herpesdrg',
        required=True,
        dest='herp',
        help='Path to herpesdrg-db.tsv'
    )
    parser.add_argument(
        '--genome', '-g',
        required=True,
        help='Path to AD169 FASTA genome'
    )
    parser.add_argument(
        '--genbank', '-b',
        required=True,
        help='Path to AD169 GenBank annotation'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output BED file path'
    )

    args = parser.parse_args()

    target_genes = ["UL97", "UL54", "UL27", "UL56", "UL51", "UL89"]
    mapped, unmapped = build_database(
        args.herp,
        args.genome,
        args.genbank,
        args.output,
        target_genes,
    )


if __name__ == "__main__":
    main()
