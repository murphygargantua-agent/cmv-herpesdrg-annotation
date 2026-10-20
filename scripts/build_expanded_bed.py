#!/usr/bin/env python3
"""Build a 3bp-spanning BED from a 1bp HerpesDRG BED for e2e variant calling."""
import sys
from pathlib import Path
import subprocess

INPUT = "data/herpesdrg_cmv_refalt.bed"
OUTPUT = "e2e_benchmark.bed"

def main():
    input_path = sys.argv[1] if len(sys.argv) > 1 else INPUT
    output_path = sys.argv[2] if len(sys.argv) > 2 else OUTPUT
    with open(input_path) as f:
        header = [l for l in f if l.startswith('#')]
        data = [l for l in f if not l.startswith('#')]
    out_gz = Path(output_path).with_suffix('.gz')
    tbi = out_gz.with_suffix('.tbi')
    for p in [out_gz, tbi]:
        if p.exists():
            p.unlink(missing_ok=True)
    with open(output_path, 'w') as f:
        for l in header:
            f.write(l)
        for l in data:
            parts = l.strip().split('\t')
            if len(parts) < 6:
                continue
            chrom, start, end, name, ref, alt = parts
            end = max(int(end), int(start) + 3)
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{ref}\t{alt}\n")
    subprocess.run(["bgzip", "-f", output_path], capture_output=True, check=True)
    subprocess.run(["tabix", "-p", "bed", str(out_gz)], capture_output=True, check=True)
    print(f"Created {output_path}.gz")

if __name__ == '__main__':
    main()
