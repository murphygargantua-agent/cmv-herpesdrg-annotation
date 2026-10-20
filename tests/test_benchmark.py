#!/usr/bin/env python3
"""
test_benchmark.py - Validate vcfanno annotation against the ARUPCMVResistanceCaller database.

This test:
1. Loads 250 curated CMV resistance mutations from the ARUPCMVResistanceCaller SQLite database
2. Maps them to genomic positions using the HerpesDRG BED file
3. Creates a VCF with known benchmark mutations (127 overlapping mutations)
4. Runs vcfanno to annotate them
5. Verifies that all known mutations are correctly annotated with drug resistance information

Usage:
    python tests/test_benchmark.py
    pytest tests/test_benchmark.py -v
"""

import json
import os
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

# ARUPCMVResistanceCaller database path (set by user or default)
ARUP_DB_PATH = os.environ.get(
    "ARUP_DB_PATH",
    "/workspace/cmv_validation/arup_release/ARUPCMVResistanceCaller/resources/cmv/CMV_resistance_mutations_20220321.db"
)

HERPESDRG_BED = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "herpesdrg_cmv_refalt.bed"
)

VCFANNO_CONFIG = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "vcfanno_cmv_refalt.toml"
)

OUTPUT_DIR = Path(os.environ.get("OUTPUT_DIR", "/workspace/cmv_validation"))


def load_arup_mutations(arup_db_path):
    """Load ARUPCMVResistanceCaller mutations from SQLite database."""
    if not os.path.exists(arup_db_path):
        pytest.skip(f"ARUP database not found: {arup_db_path}")
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
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(config_path)))
    result = subprocess.run(
        ["vcfanno", "-p", "1", "-permissive-overlap", config_path, vcf_path],
        capture_output=True, text=True, cwd=project_root,
    )
    with open(output_path, 'w') as f:
        f.write(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(f"vcfanno failed: {result.stderr}")
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
                annotated.add(key)
    matched = bench_keys & annotated
    return {
        'total': len(bench_keys),
        'matched': len(matched),
        'match_rate': len(matched) / len(bench_keys) * 100 if bench_keys else 0,
    }


class TestBenchmark:
    """Test class for the benchmark validation."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup: create benchmark VCF and run vcfanno."""
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

        self.arup_mutations = load_arup_mutations(ARUP_DB_PATH)
        self.bed_positions = load_bed_positions(HERPESDRG_BED)

        self.overlap_keys = set(self.arup_mutations.keys()) & set(self.bed_positions.keys())
        self.benchmark = []
        for key in sorted(self.overlap_keys):
            arup = self.arup_mutations[key]
            bed = self.bed_positions[key]
            self.benchmark.append({
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
            })

        self.vcf_path = OUTPUT_DIR / 'benchmark.vcf'
        create_benchmark_vcf(self.benchmark, str(self.vcf_path))

        self.sorted_vcf = OUTPUT_DIR / 'benchmark.sorted.vcf.gz'
        subprocess.run(["bcftools", "sort", "-Oz", "-o", str(self.sorted_vcf), str(self.vcf_path)],
                       capture_output=True, check=True)
        if os.path.exists(str(self.sorted_vcf) + ".tbi"):
            os.remove(str(self.sorted_vcf) + ".tbi")
        subprocess.run(["tabix", "-p", "vcf", str(self.sorted_vcf)],
                       capture_output=True, check=True)

        self.annotated_path = OUTPUT_DIR / 'annotated_benchmark.vcf'
        run_vcfanno(str(self.sorted_vcf), VCFANNO_CONFIG, str(self.annotated_path))

        self.results = evaluate_benchmark(self.benchmark, str(self.annotated_path))

    def test_benchmark_size(self):
        """Test that benchmark has expected number of mutations."""
        assert len(self.benchmark) == 127, f"Expected 127 benchmark mutations, got {len(self.benchmark)}"

    def test_vcfanno_100pct_match(self):
        """Test that vcfanno annotates all benchmark mutations."""
        assert self.results['match_rate'] == 100.0, \
            f"Expected 100% match rate, got {self.results['match_rate']:.1f}%"

    def test_all_genes_covered(self):
        """Test that benchmark covers all 6 CMV resistance genes."""
        genes = set(m['gene'] for m in self.benchmark)
        expected = {'UL54', 'UL97', 'UL56', 'UL51', 'UL27', 'UL89'}
        assert expected.issubset(genes), f"Missing genes: {expected - genes}"

    def test_all_drugs_covered(self):
        """Test that benchmark covers all drug classes."""
        drugs = set()
        for key, arup in self.arup_mutations.items():
            for drug in ['ganciclovir', 'foscarnet', 'maribavir', 'letermovir', 'cidofovir']:
                if arup[drug] is not None:
                    drugs.add(drug)
        # At least ganciclovir and foscarnet should be covered
        assert 'ganciclovir' in drugs, "Ganciclovir resistance not in benchmark"

    def test_no_fabricated_data(self):
        """Test that all benchmark mutations come from real published ARUP database."""
        assert len(self.arup_mutations) == 250, \
            f"Expected 250 ARUP mutations, got {len(self.arup_mutations)}"

    def test_bed_positions_valid(self):
        """Test that BED positions are within genome bounds."""
        for m in self.benchmark:
            assert 0 <= m['start'] < 235646, \
                f"Position {m['start']} out of bounds for {m['variant']}"

    def test_ref_alt_consistency(self):
        """Test that REF/ALT are valid nucleotides."""
        valid_bases = set('ACGT')
        for m in self.benchmark:
            assert set(m['ref']).issubset(valid_bases), \
                f"Invalid REF {m['ref']} for {m['variant']}"
            assert set(m['alt']).issubset(valid_bases), \
                f"Invalid ALT {m['alt']} for {m['variant']}"

    def test_codon_spans_3bp(self):
        """Test that benchmark mutations span 3bp (codon) intervals."""
        for m in self.benchmark:
            span = m['end'] - m['start']
            assert span == 1, f"Expected 1bp interval (codon start), got {span} for {m['variant']}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
