#!/usr/bin/env python3
"""Tests for HerpesDRG vcfanno database."""
import os
import pysam


def test_bed_gz_exists():
    path = "data/herpesdrg_cmv_refalt.bed.gz"
    assert os.path.exists(path), f"{path} does not exist"


def test_tabix_index_exists():
    path = "data/herpesdrg_cmv_refalt.bed.gz.tbi"
    assert os.path.exists(path), f"{path} does not exist"


def test_chromosome_name():
    tb = pysam.TabixFile("data/herpesdrg_cmv_refalt.bed.gz")
    assert "NC_006273.2" in tb.contigs, f"Expected NC_006273.2 in contigs, got {tb.contigs}"


def test_total_entries():
    tb = pysam.TabixFile("data/herpesdrg_cmv_refalt.bed.gz")
    count = sum(1 for _ in tb.fetch("NC_006273.2", 0, 235646))
    assert count > 300, f"Expected >300 entries, got {count}"


def test_has_ul54_d301n():
    """Test UL54_D301N with correct REF/ALT."""
    tb = pysam.TabixFile("data/herpesdrg_cmv_refalt.bed.gz")
    found = False
    for line in tb.fetch("NC_006273.2", 81019, 81022):
        if "UL54_D301N" in line:
            found = True
            parts = line.rstrip().split("\t")
            assert parts[4] == "ATC", f"Expected REF=ATC, got {parts[4]}"
            assert parts[5] == "ATT", f"Expected ALT=ATT, got {parts[5]}"
    assert found, "UL54_D301N not found"


def test_has_ul27_r448p():
    """Test UL27_R448P (a known HerpesDRG mutation)."""
    tb = pysam.TabixFile("data/herpesdrg_cmv_refalt.bed.gz")
    found = False
    for line in tb.fetch("NC_006273.2", 33663, 33666):
        if "UL27_R448P" in line:
            found = True
            parts = line.rstrip().split("\t")
            # R448P is at position 33664-33665 in the original BED
            # Verify it has 6 columns
            assert len(parts) == 6, f"Expected 6 columns, got {len(parts)}"
    assert found, "UL27_R448P not found"


def test_refalt_header():
    with open("data/herpesdrg_cmv_refalt.bed", "r") as f:
        header = f.readline()
    assert "chr" in header, "Header should contain chr"
    assert "ref" in header, "Header should contain ref"
    assert "alt" in header, "Header should contain alt"


def test_refalt_columns():
    with open("data/herpesdrg_cmv_refalt.bed", "r") as f:
        f.readline()
        for line in f:
            parts = line.rstrip().split("\t")
            if not line.startswith("#"):
                assert len(parts) == 6, f"Expected 6 columns, got {len(parts)}"


def test_positions_in_nc_006273():
    tb = pysam.TabixFile("data/herpesdrg_cmv_refalt.bed.gz")
    max_pos = 0
    for line in tb.fetch("NC_006273.2", 0, 235646):
        start = int(line.split("\t")[1])
        max_pos = max(max_pos, start)
    assert max_pos < 235646, f"Position {max_pos} exceeds genome length"


def test():
    tests = [
        test_bed_gz_exists,
        test_tabix_index_exists,
        test_chromosome_name,
        test_total_entries,
        test_has_ul54_d301n,
        test_has_ul27_r448p,
        test_refalt_header,
        test_refalt_columns,
        test_positions_in_nc_006273,
    ]
    for test in tests:
        try:
            test()
            print(f"  PASS: {test.__name__}")
        except Exception as e:
            print(f"  FAIL: {test.__name__}: {e}")


if __name__ == "__main__":
    test()
