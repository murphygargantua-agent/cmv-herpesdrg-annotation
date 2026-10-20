#!/usr/bin/env python3
"""
End-to-end (e2e) benchmark test for the HerpesDRG vcfanno database.

Simulates a real CMV NGS workflow:
1. Takes known mutations from the ARUPCMVResistanceCaller database
2. Creates synthetic FASTQ reads with high coverage at mutation positions
3. Aligns reads with BWA
4. Calls variants with bcftools
5. Annotates variants with the HerpesDRG vcfanno database (3bp-spanning BED)
6. Verifies that all known mutations are correctly detected and annotated
"""

import json
import os
import random
import sqlite3
import subprocess
from pathlib import Path

import pysam
import pytest

ARUP_DB_PATH = os.environ.get(
    "ARUP_DB_PATH",
    "/workspace/cmv_validation/arup_release/ARUPCMVResistanceCaller/resources/cmv/CMV_resistance_mutations_20220321.db"
)
REF_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "AD169_NC_006273.2.fa"
)
BED_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "herpesdrg_cmv_refalt.bed.gz"
)
EXPANDED_BED = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "e2e_benchmark.bed.gz")
OUTPUT_DIR = Path(os.environ.get("OUTPUT_DIR", "/workspace/cmv_validation/e2e_test"))


def load_arup_mutations():
    if not os.path.exists(ARUP_DB_PATH):
        pytest.skip(f"ARUP database not found: {ARUP_DB_PATH}")
    db = sqlite3.connect(ARUP_DB_PATH)
    c = db.cursor()
    c.execute("SELECT gene, variant, codon_start, codon_end, ref_aa, var_aa, "
              "cidofovir, foscarnet, ganciclovir, maribavir, letermovir FROM mutations")
    mutations = []
    for row in c.fetchall():
        mutations.append({
            'gene': row[0], 'variant': row[1],
            'codon_start': int(row[2]), 'codon_end': int(row[3]),
            'ref_aa': row[4], 'var_aa': row[5],
            'ganciclovir': row[6], 'foscarnet': row[7],
            'maribavir': row[8], 'letermovir': row[9], 'cidofovir': row[10],
        })
    db.close()
    return mutations


def load_bed_positions():
    positions = {}
    with open(BED_PATH.replace('.gz', '')) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            chrom, start, end, name, ref, alt = parts
            base = name.split('|')[0]
            gene, variant = base.split('_', 1)
            positions[f"{gene}|{variant}"] = {
                'chrom': chrom, 'start': int(start), 'end': int(end),
                'ref': ref, 'alt': alt, 'gene': gene, 'variant': variant
            }
    return positions


def create_synthetic_reads(mutations, ref_path, n_reads=5000, read_len=150):
    with open(ref_path) as f:
        ref_name = f.readline().strip()[1:]
        ref_seq = ''.join(line.strip() for line in f)
    ref_list = list(ref_seq)

    bed_positions = load_bed_positions()
    selected = []
    for m in mutations:
        key = f"{m['gene']}|{m['variant']}"
        if key in bed_positions:
            selected.append(m)
            bed = bed_positions[key]
            for i, base in enumerate(bed['alt']):
                pos = bed['start'] + i
                if pos < len(ref_list):
                    ref_list[pos] = base

    if not selected:
        pytest.skip("No mutations found in BED overlap")

    ref_str_mutated = ''.join(ref_list)
    reads = []

    # 50 reads covering each mutation position
    for m in selected:
        key = f"{m['gene']}|{m['variant']}"
        bed = bed_positions[key]
        for _ in range(50):
            start = max(0, bed['start'] - 75)
            read = ref_str_mutated[start:start + read_len]
            if len(read) == read_len:
                reads.append(read)

    # Background reads
    random.seed(42)
    n_background = n_reads - len(reads)
    for _ in range(n_background):
        start = random.randint(0, max(0, len(ref_str_mutated) - read_len))
        read = ref_str_mutated[start:start + read_len]
        reads.append(read)

    random.shuffle(reads)
    fastq_path = OUTPUT_DIR / "e2e_test.fastq"
    with open(fastq_path, 'w') as f:
        for i, read in enumerate(reads):
            f.write(f"@read_{i}\n{read}\n+\n{'I' * len(read)}\n")
    return selected, fastq_path


def run_alignment(fastq_path, ref_path, bam_path):
    subprocess.run(["bwa", "index", ref_path], capture_output=True, check=True)
    # Ensure expanded BED is indexed
    expanded_bed = EXPANDED_BED.replace('.gz', '')
    if not os.path.exists(expanded_bed + '.gz'):
        subprocess.run(["bgzip", "-f", expanded_bed], capture_output=True, check=True)
    subprocess.run(["tabix", "-p", "bed", expanded_bed + '.gz'],
                   capture_output=True, check=True)
    sam_path = OUTPUT_DIR / "e2e_test.sam"
    with open(str(sam_path), 'w') as f:
        subprocess.run(["bwa", "mem", "-t", "4", ref_path, fastq_path],
                       check=True, stdout=f, stderr=subprocess.DEVNULL)
    subprocess.run(["samtools", "sort", "-o", str(bam_path), str(sam_path)],
                   capture_output=True, check=True)
    subprocess.run(["samtools", "index", str(bam_path)], capture_output=True, check=True)
    if sam_path.exists():
        os.remove(sam_path)


def run_variant_calling(bam_path, vcf_path):
    mp = subprocess.run(
        ["bcftools", "mpileup", "-f", REF_PATH, "--max-depth", "10000",
         "-Ou", str(bam_path)],
        capture_output=True, check=True,
    )
    subprocess.run(
        ["bcftools", "call", "-mv", "-Oz", "-o", str(vcf_path)],
        input=mp.stdout, check=True,
    )
    subprocess.run(["tabix", "-p", "vcf", str(vcf_path)],
                   capture_output=True, check=True)


def run_vcfanno(vcf_path, output_path):
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    config_path = OUTPUT_DIR / "e2e_config.toml"
    with open(str(config_path), 'w') as f:
        f.write('[[annotation]]\n')
        f.write(f'file = "{EXPANDED_BED}"\n')
        f.write('columns = [4]\n')
        f.write('names = ["cmv_resistance"]\n')
        f.write('ops = ["uniq"]\n')
        f.write('type = "String"\n')
    result = subprocess.run(
        ["vcfanno", "-p", "1", "-permissive-overlap", str(config_path), str(vcf_path)],
        capture_output=True, text=True, cwd=project_root,
    )
    with open(str(output_path), 'w') as f:
        f.write(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(f"vcfanno failed: {result.stderr[:500]}")


def evaluate_e2e(benchmark_mutations, annotated_path):
    bench_variants = set(m['variant'] for m in benchmark_mutations)
    detected = set()
    with open(str(annotated_path)) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            info = parts[7]
            if 'cmv_resistance=' in info:
                ann = info.split('cmv_resistance=')[1].split(';')[0]
                mutation_name = ann.split('|')[0]
                gene, variant = mutation_name.split('_', 1)
                detected.add(variant)
    matched = bench_variants & detected
    return {
        'total': len(bench_variants),
        'detected': len(detected & bench_variants),
        'matched': len(matched),
        'match_rate': len(matched) / len(bench_variants) * 100 if bench_variants else 0,
        'bench_variants': sorted(bench_variants),
        'detected_variants': sorted(detected & bench_variants),
    }


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    arup_mutations = load_arup_mutations()
    print(f"ARUP mutations: {len(arup_mutations)}")

    bed_positions = load_bed_positions()
    gene_mutations = {}
    for m in arup_mutations:
        key = f"{m['gene']}|{m['variant']}"
        if key in bed_positions:
            gene_mutations.setdefault(m['gene'], []).append(m)

    selected = []
    for gene in sorted(gene_mutations.keys()):
        selected.extend(gene_mutations[gene][:5])
    print(f"Selected {len(selected)} benchmark mutations")

    n_reads = 5000
    selected, fastq_path = create_synthetic_reads(selected, REF_PATH, n_reads)
    bam_path = OUTPUT_DIR / "e2e_test.bam"
    vcf_path = OUTPUT_DIR / "e2e_test.vcf.gz"
    annotated_path = OUTPUT_DIR / "e2e_annotated.vcf"

    print(f"Creating synthetic FASTQ with {n_reads} reads...")
    print("Aligning with BWA...")
    run_alignment(fastq_path, REF_PATH, bam_path)
    print("Calling variants with bcftools...")
    run_variant_calling(bam_path, vcf_path)
    print("Annotating with HerpesDRG vcfanno...")
    run_vcfanno(vcf_path, annotated_path)

    results = evaluate_e2e(selected, annotated_path)
    print(f"\n=== E2E Results ===")
    print(f"Benchmark mutations: {results['total']}")
    print(f"Detected: {results['detected']}")
    print(f"Match rate: {results['match_rate']:.1f}%")
    if results['match_rate'] == 100.0:
        print("All benchmark mutations correctly detected and annotated!")
    else:
        print(f"{results['total'] - results['detected']} mutations not detected")

    results_path = OUTPUT_DIR / "e2e_benchmark_results.json"
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to: {results_path}")


if __name__ == '__main__':
    main()
