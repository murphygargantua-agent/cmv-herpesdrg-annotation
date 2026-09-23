#!/usr/bin/env python3
"""
build_database.py - Rebuild the CMV HerpesDRG vcfanno database.

Builds codon-aware BED entries by mapping amino-acid mutations to their
genomic codon coordinates in the AD169 reference genome, using the
aa_to_genomic mapping derived from the GenBank CDS translations.

Usage:
    python build_database.py --herpesdrg herpesdrg-db.tsv \
                             --genome AD169.fasta \
                             --genbank AD169_annotation.gb \
                             --output herpesdrg_cmv.bed

Downloads:
    - HerpesDRG: https://github.com/ojcharles/herpesdrg-db/raw/main/herpesdrg-db.tsv
    - AD169 genome: NCBI accession X17403.1 (via Entrez EFetch)
    - AD169 annotation: NCBI accession X17403.1 (via Entrez EFetch)
"""

import argparse
import csv
import json
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq


def parse_fasta(path: str) -> str:
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                seq.append(line)
    return ''.join(seq)


def parse_genbank_cds(target_genes: list, genbank_path: str) -> dict:
    """Parse GenBank CDS for target genes. Returns dict with gene name -> info."""
    with open(genbank_path) as f:
        record = SeqIO.read(f, 'genbank')

    gene_info = {}
    for feature in record.features:
        if feature.type == "CDS" and "gene" in feature.qualifiers:
            raw_gene = feature.qualifiers["gene"][0]
            gene_name = raw_gene[4:].split()[0] if raw_gene.startswith("HCMV") else raw_gene.split()[0]
            if gene_name in target_genes:
                start0 = int(feature.location.start)  # 0-based
                end0 = int(feature.location.end)      # 0-based exclusive
                strand = feature.location.strand
                gene_info[gene_name] = {
                    'start0': start0,
                    'end0': end0,
                    'strand': strand,
                    'length': end0 - start0,
                }
    return gene_info


def codon_start_0based(gene_info: dict, aa_pos: int) -> int:
    """
    Compute 0-based genomic start of the codon for amino acid position P (1-based).

    For forward strand: codon = gene_start0 + 3*(P-1)
    For reverse strand: codon = gene_end0 - 3*P
    """
    start0 = gene_info['start0']
    end0 = gene_info['end0']
    strand = gene_info['strand']

    if strand == 1:
        return start0 + 3 * (aa_pos - 1)
    else:
        return end0 - 3 * aa_pos


def verify_codon(genome: str, gene_info: dict, aa_pos: int, ref_aa: str) -> bool:
    """Extract codon and verify it translates to ref_aa."""
    start0 = codon_start_0based(gene_info, aa_pos)
    if start0 < 0 or start0 + 3 > len(genome):
        return False
    codon_plus = genome[start0:start0 + 3]
    if gene_info['strand'] == -1:
        codon = str(Seq(codon_plus).reverse_complement())
    else:
        codon = codon_plus
    return str(Seq(codon).translate()) == ref_aa


def build_database(herpesdrg_path: str, genome_path: str,
                   genbank_path: str, output_path: str,
                   target_genes: list) -> tuple:
    genome = parse_fasta(genome_path)
    print(f"Genome: {len(genome)} bp")

    gene_info = parse_genbank_cds(target_genes, genbank_path)
    print(f"Found {len(gene_info)} genes:")
    for g in sorted(gene_info.keys()):
        info = gene_info[g]
        print(f"  {g}: 0-based {info['start0']}-{info['end0']} ({info['length']} bp) strand={info['strand']}")

    # Load aa_to_genomic mapping
    mapping_path = os.path.join(os.path.dirname(genbank_path), "ad169_gene_mapping.json")
    if not os.path.exists(mapping_path):
        mapping_path = "data/ad169_gene_mapping.json"
    if not os.path.exists(mapping_path):
        raise FileNotFoundError(f"aa_to_genomic mapping not found at {mapping_path}")
    with open(mapping_path) as f:
        aa_map = json.load(f)

    print("\nLoading HerpesDRG TSV...")
    with open(herpesdrg_path) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    print(f"Total rows: {len(rows)}")
    target_rows = [r for r in rows if r.get('gene') in target_genes]
    print(f"Target gene rows: {len(target_rows)}")

    print("\nBuilding mutation list...")
    mutations = []
    unmapped = 0
    unparsed = 0
    verified = 0
    failed_verify = 0

    for row in target_rows:
        gene = row['gene']
        aa_change = row['aa_change']

        match = re.match(r'^([A-Z])(\d+)([A-Z])$', aa_change)
        if not match:
            unparsed += 1
            continue
        ref_aa, pos, mut_aa = match.group(1), int(match.group(2)), match.group(3)

        if gene not in gene_info:
            unmapped += 1
            continue

        if pos < 1 or pos > gene_info[gene]['length']:
            unmapped += 1
            continue

        if not verify_codon(genome, gene_info[gene], pos, ref_aa):
            failed_verify += 1
            continue

        start0 = codon_start_0based(gene_info[gene], pos)
        chrom = 'NC_006273'

        drug_info = []
        for drug in ['Ganciclovir', 'Aciclovir', 'Cidofovir', 'Foscarnet',
                     'Brincidofovir', 'Letermovir', 'Maribavir']:
            val = row.get(drug, '')
            if val and str(val).strip():
                drug_info.append(f"{drug}={val}")

        # Label: GENE_AA_CHANGE|Drug1=val1;Drug2=val2
        label = f"{gene}_{aa_change}"
        if drug_info:
            label += '|' + ';'.join(drug_info)

        # Bed entry: codon-spanning 3-base interval (start0, start0+3)
        mutations.append({
            'chrom': chrom,
            'start': start0,
            'end': start0 + 3,
            'label': label,
        })
        verified += 1

    # Sort by chrom then start
    mutations.sort(key=lambda m: (m['chrom'], m['start']))

    # Merge entries at the same genomic position (same codon, different mutations/drugs)
    merged = {}
    for m in mutations:
        key = (m['chrom'], m['start'], m['end'])
        if key not in merged:
            merged[key] = [m['label']]
        else:
            merged[key].append(m['label'])
    # Deduplicate within each key
    for key in merged:
        seen = set()
        labels = []
        for l in merged[key]:
            if l not in seen:
                seen.add(l)
                labels.append(l)
        merged[key] = labels

    print(f"\nMapped/verified: {verified}")
    print(f"Unmapped: {unmapped}")
    print(f"Unparsed: {unparsed}")
    print(f"Failed codon verify: {failed_verify}")
    print(f"Unique positions: {len(merged)} (from {len(mutations)} entries)")

    with open(output_path, 'w') as f:
        f.write("##bedFormat=4\n")
        for (chrom, start, end), labels in sorted(merged.items(), key=lambda x: (x[0][0], x[0][1])):
            f.write(f"{chrom}\t{start}\t{end}\t{';'.join(labels)}\n")

    print(f"\nWrote {len(merged)} unique positions to {output_path}")
    return len(merged), unmapped + unparsed + failed_verify


def main():
    parser = argparse.ArgumentParser(description="Build CMV HerpesDRG vcfanno database")
    parser.add_argument("--herpesdrg", required=True, help="Path to herpesdrg-db.tsv")
    parser.add_argument("--genome", required=True, help="Path to AD169 FASTA")
    parser.add_argument("--genbank", required=True, help="Path to AD169 GenBank annotation")
    parser.add_argument("--output", default="herpesdrg_cmv.bed", help="Output BED file")
    parser.add_argument("--genes", nargs='+', default=["UL54","UL97","UL56","UL51","UL27","UL89"])
    args = parser.parse_args()

    build_database(args.herpesdrg, args.genome, args.genbank, args.output, args.genes)


if __name__ == "__main__":
    main()
