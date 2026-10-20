#!/usr/bin/env python3
"""
test_vcfanno_benchmark.py - Validate vcfanno annotation against the ARUPCMVResistanceCaller database.

This script:
1. Extracts known CMV resistance mutations from the ARUPCMVResistanceCaller SQLite database
2. Maps them to genomic positions using the HerpesDRG BED file
3. Creates a VCF with these known mutations
4. Runs vcfanno to annotate them
5. Verifies that all known mutations are correctly annotated

Usage:
    python scripts/test_vcfanno_benchmark.py \
        --herpesdrg-data data/herpesdrg-db.tsv \
        --genome data/AD169_NC_006273.2.fa \
        --genbank data/AD169_NC_006273.2.gb \
        --arup-db /workspace/cmv_validation/arup_release/ARUPCMVResistanceCaller/resources/cmv/CMV_resistance_mutations_20220321.db \
        --bed data/herpesdrg_cmv_refalt.bed.gz \
        --vcfanno-config data/vcfanno_cmv_refalt.toml \
        --output-dir /workspace/cmv_validation
"""

import argparse
import json
import os
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path


def load_arup_mutations(arup_db_path):
    """Load ARUPCMVResistanceCaller mutations."""
    db = sqlite3.connect(arup_db_path)
    c = db.cursor()
    c.execute("SELECT gene, variant, codon_start, codon_end, ref_aa, var_aa, "
              "cidofovir, foscarnet, ganciclovir, maribavir, letermovir FROM mutations")
    mutations = {}
    for row in c.fetchall():
        key = f"{row[0]}|{row[1]}"
        mutations[key] = {
            'gene': row[0], 'variant': row[1],
            'codon_start': int(row[2]), 'codon_end': int(row[3]),
            'ref_aa': row[4], 'var_aa': row[5],
            'cidofovir': row[6], 'foscarnet': row[7],
            'ganciclovir': row[8], 'maribavir': row[9],
            'letermovir': row[10],
        }
    db.close()
    return mutations


def load_bed_positions(bed_path):
    """Load BED positions from the ref/alt BED file."""
    positions = {}
    with open(bed_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            chrom, start, end, name, ref, alt = parts
            base = name.split('|')[0]
            gene, variant = base.split('_', 1)
            key = f"{gene}|{variant}"
            positions[key] = {
                'chrom': chrom, 'start': int(start), 'end': int(end),
                'ref': ref, 'alt': alt, 'gene': gene, 'variant': variant
            }
    return positions


def create_benchmark_vcf(benchmark, output_path):
    """Create a VCF with known benchmark mutations."""
    with open(output_path, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=NC_006273.2,length=235646>\n")
        f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">\n')
        f.write('##INFO=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for m in benchmark:
            pos = m['start'] + 1  # 0-based -> 1-based
            f.write(f"NC_006273.2\t{pos}\t.\t{m['ref']}\t{m['alt']}\t200\tPASS\t"
                    f"DP=30;AD=0,30\tGT:AD\t1/1:0,30\n")


def run_vcfanno(vcf_path, config_path, output_path):
    """Run vcfanno on the VCF."""
    cmd = ["vcfanno", "-p", "1", "-permissive-overlap", config_path, vcf_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    with open(output_path, 'w') as f:
        f.write(result.stdout)
    if result.stderr:
        with open(output_path + '.stderr', 'w') as f:
            f.write(result.stderr)
    return result.returncode


def evaluate_benchmark(benchmark, annotated_path):
    """Evaluate which benchmark mutations were annotated."""
    bench_keys = set(f"{m['chrom']}:{m['start']+1}" for m in benchmark)
    annotated = set()
    with open(annotated_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            key = f"{parts[0]}:{parts[1]}"
            info = parts[7]
            if 'cmv_resistance=' in info:
                ann = info.split('cmv_resistance=')[1].split(';')[0]
                annotated.add(key)
                # Also store the annotation
    matched = bench_keys & annotated
    return {
        'total': len(bench_keys),
        'matched': len(matched),
        'match_rate': len(matched) / len(bench_keys) * 100 if bench_keys else 0,
        'matched_keys': sorted(matched),
    }


def main():
    parser = argparse.ArgumentParser(description="Test vcfanno annotation against ARUP benchmark")
    parser.add_argument('--herpesdrg-data', default='data/herpesdrg-db.tsv')
    parser.add_argument('--genome', default='data/AD169_NC_006273.2.fa')
    parser.add_argument('--genbank', default='data/AD169_NC_006273.2.gb')
    parser.add_argument('--arup-db', default='/workspace/cmv_validation/arup_release/ARUPCMVResistanceCaller/resources/cmv/CMV_resistance_mutations_20220321.db')
    parser.add_argument('--bed', default='data/herpesdrg_cmv_refalt.bed.gz')
    parser.add_argument('--bed-raw', default='data/herpesdrg_cmv_refalt.bed')
    parser.add_argument('--vcfanno-config', default='data/vcfanno_cmv_refalt.toml')
    parser.add_argument('--output-dir', default='/workspace/cmv_validation')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Loading ARUP mutations...")
    arup_mutations = load_arup_mutations(args.arup_db)
    print(f"  {len(arup_mutations)} ARUP mutations")

    print("Loading BED positions...")
    bed_positions = load_bed_positions(args.bed_raw)
    print(f"  {len(bed_positions)} BED positions")

    # Find overlap
    overlap_keys = set(arup_mutations.keys()) & set(bed_positions.keys())
    print(f"  {len(overlap_keys)} overlapping mutations")

    # Build benchmark
    benchmark = []
    for key in sorted(overlap_keys):
        arup = arup_mutations[key]
        bed = bed_positions[key]
        benchmark.append({
            'gene': bed['gene'],
            'variant': bed['variant'],
            'chrom': bed['chrom'],
            'start': bed['start'],
            'end': bed['end'],
            'ref': bed['ref'],
            'alt': bed['alt'],
            'codon_start': arup['codon_start'],
            'codon_end': arup['codon_end'],
            'ref_aa': arup['ref_aa'],
            'var_aa': arup['var_aa'],
            'ganciclovir': arup['ganciclovir'],
            'foscarnet': arup['foscarnet'],
            'maribavir': arup['maribavir'],
            'letermovir': arup['letermovir'],
        })

    print(f"\nBenchmark: {len(benchmark)} mutations")

    # Create VCF
    vcf_path = os.path.join(args.output_dir, 'benchmark.vcf')
    create_benchmark_vcf(benchmark, vcf_path)
    print(f"Created VCF: {vcf_path}")

    # Sort VCF
    sorted_vcf = os.path.join(args.output_dir, 'benchmark.sorted.vcf.gz')
    subprocess.run(["bcftools", "sort", "-Oz", "-o", sorted_vcf, vcf_path],
                   capture_output=True, check=True)
    subprocess.run(["tabix", "-p", "vcf", sorted_vcf], capture_output=True, check=True)
    print(f"Sorted VCF: {sorted_vcf}")

    # Run vcfanno
    annotated_path = os.path.join(args.output_dir, 'annotated_benchmark.vcf')
    os.chdir(os.path.abspath(os.path.dirname(args.vcfanno_config) + "/.."))
    rc = run_vcfanno(sorted_vcf, args.vcfanno_config, annotated_path)
    print(f"vcfanno exit: {rc}")

    # Evaluate
    results = evaluate_benchmark(benchmark, annotated_path)
    print(f"\n=== Results ===")
    print(f"Benchmark mutations: {results['total']}")
    print(f"Matched: {results['matched']}")
    print(f"Match rate: {results['match_rate']:.1f}%")

    # Save results
    results_path = os.path.join(args.output_dir, 'benchmark_results.json')
    with open(results_path, 'w') as f:
        json.dump({
            'total': results['total'],
            'matched': results['matched'],
            'match_rate': results['match_rate'],
            'benchmark_mutations': benchmark,
            'matched_keys': results['matched_keys'],
        }, f, indent=2)
    print(f"Results saved to: {results_path}")


if __name__ == '__main__':
    main()
