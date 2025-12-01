#!/usr/bin/env python3

"""
Group HHR files by overlapping coordinates and strand.
"""

import os, re, sys, argparse
from collections import defaultdict, Counter
from pathlib import Path

# Regex: start-end_(+/-)_frame.hhr

PAT = re.compile(r'^(\d+)-(\d+)_\(([+-])\)_(\d+)\.hhr$')

def find_files(hhrdir):

    """Find all .hhr files in directory."""

    files = []
    for name in os.listdir(hhrdir):
        if name.endswith(".hhr"):
            files.append(os.path.join(hhrdir, name))
    return sorted(files)

class DSU:

    """Union-Find / Disjoint Set Union data structure."""

    def __init__(self, n):
        self.p = list(range(n))
        self.r = [0] * n
    
    def f(self, x):

        """Find with path compression."""

        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x
    
    def u(self, a, b):

        """Union by rank."""

        ra, rb = self.f(a), self.f(b)
        if ra == rb:
            return
        if self.r[ra] < self.r[rb]:
            ra, rb = rb, ra
        self.p[rb] = ra
        if self.r[ra] == self.r[rb]:
            self.r[ra] += 1

def parse_args():

    """Parse command line arguments."""

    parser = argparse.ArgumentParser(
        description='Group HHR files by overlapping coordinates and strand'
    )
    parser.add_argument(
        '--hhrdir',
        required=True,
        help='Directory containing .hhr files'
    )
    parser.add_argument(
        '--outdir',
        default='.',
        help='Output directory (default: current directory)'
    )
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Validate input directory

    if not os.path.isdir(args.hhrdir):
        print(f"ERROR: Directory not found: {args.hhrdir}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory

    os.makedirs(args.outdir, exist_ok=True)
    
    # Find files

    files = find_files(args.hhrdir)
    recs = []  # (path, start, end, strand, frame)
    skipped = []

    for path in files:
        base = os.path.basename(path)
        m = PAT.match(base)
        if not m:
            skipped.append(base)
            continue
        start = int(m.group(1))
        end = int(m.group(2))
        strand = m.group(3)
        frame = int(m.group(4))
        recs.append((path, start, end, strand, frame))

    n = len(recs)
    if n == 0:
        print("WARNING: No .hhr files with expected format found.", file=sys.stderr)
        print("Expected format: start-end_(+/-)_frame.hhr", file=sys.stderr)
        if skipped:
            print(f"Files with unexpected format (skipped): {len(skipped)}", file=sys.stderr)
            for s in skipped[:10]:  # Show first 10
                print(f"  {s}", file=sys.stderr)
        sys.exit(0)

    # Build DSU - connect by (strand, start) and (strand, end)

    dsu = DSU(n)
    key_start = {}  # (strand, start) -> idx
    key_end = {}    # (strand, end)   -> idx

    for i, (path, start, end, strand, frame) in enumerate(recs):
        ks = (strand, start)
        ke = (strand, end)
        if ks in key_start:
            dsu.u(i, key_start[ks])
        else:
            key_start[ks] = i
        if ke in key_end:
            dsu.u(i, key_end[ke])
        else:
            key_end[ke] = i

    # Group by component

    comps = defaultdict(list)
    for i, r in enumerate(recs):
        comps[dsu.f(i)].append(r)

    # Output files

    groups_path = os.path.join(args.outdir, "groups.tsv")
    singles_path = os.path.join(args.outdir, "singles.tsv")
    summary_path = os.path.join(args.outdir, "groups_by_size_summary.tsv")
    skipped_path = os.path.join(args.outdir, "skipped_unexpected_format.txt")

    sizes_counter = Counter()
    total_groups = 0
    total_members = 0

    with open(groups_path, "w") as g, open(singles_path, "w") as sng:
        g.write("group_id\tsize\tstrand\tstarts\tends\tframes\tfiles\n")
        sng.write("file\tstrand\tstart\tend\tframe\n")
        gid = 1
        for root, members in sorted(comps.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            size = len(members)
            strand_set = {m[3] for m in members}
            strand = strand_set.pop() if len(strand_set) == 1 else "MIXED"

            starts = sorted({m[1] for m in members})
            ends = sorted({m[2] for m in members})
            frames = sorted({m[4] for m in members})
            files_list = [
                os.path.basename(m[0]) 
                for m in sorted(members, key=lambda x: (x[1], x[2], x[4]))
            ]

            if size == 1:
                (path, start, end, sd, frame) = members[0]
                sng.write(f"{os.path.basename(path)}\t{sd}\t{start}\t{end}\t{frame}\n")
            else:
                g.write(
                    f"{gid}\t{size}\t{strand}\t"
                    f"{','.join(map(str, starts))}\t"
                    f"{','.join(map(str, ends))}\t"
                    f"{','.join(map(str, frames))}\t"
                    f"{';'.join(files_list)}\n"
                )
                gid += 1
                total_groups += 1
                total_members += size
                sizes_counter[size] += 1

    with open(summary_path, "w") as sm:
        sm.write("group_size\tn_groups\n")
        for sz in sorted(sizes_counter):
            sm.write(f"{sz}\t{sizes_counter[sz]}\n")

    if skipped:
        with open(skipped_path, "w") as sk:
            for x in skipped:
                sk.write(x + "\n")

    # Summary

    n_singles = n - total_members
    print("===== GROUPING SUMMARY =====")
    print(f"HHR files with expected format   : {n}")
    print(f"HHR files with unexpected format : {len(skipped)}")
    print(f"Groups (size >= 2)               : {total_groups}")
    for sz in sorted(sizes_counter):
        print(f"  - Groups of size {sz:2d}          : {sizes_counter[sz]}")
    print(f"Singles (size = 1)               : {n_singles}")
    print("\nOutput files:")
    print(f"  {groups_path}")
    print(f"  {singles_path}")
    print(f"  {summary_path}")
    if skipped:
        print(f"  {skipped_path}")

if __name__ == "__main__":
    main()

