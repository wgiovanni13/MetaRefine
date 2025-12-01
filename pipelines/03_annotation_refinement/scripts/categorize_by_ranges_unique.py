#!/usr/bin/env python3

"""
Categorize unique genomic ranges from a combined overlap table.

Exact logic preserved from the original script:
- No filtering by Feature_* (all rows are considered, as in the original).
- Same promotion rules from Cat3 -> Cat2 based on tools count and % overlap.
- Same classification priority and output file names/columns.

Optional CLI flags are provided for Nextflow integration, but defaults replicate
the original behavior and results.
"""

from __future__ import annotations
import argparse
import csv
import os
from typing import Dict, List, Set, Tuple

# ----- Defaults matching the original script -----

DEFAULT_INPUT = "/users/wagner.guzman/scripts/loki_b36/frameshift/final_overlap_results_combined.tsv"
DEFAULT_THRESH = 10.0

CAT1_FILE = "cat1_single_tool.tsv"
CAT2_FILE = "cat2_full_consensus_5.tsv"
CAT3_FILE = "cat3_partial_overlaps.tsv"
CAT4_FILE = "cat4_other_full_or_uncategorized.tsv"
SUMMARY_FILE = "category_summary.tsv"
ALL_UNIQUE_FILE = "all_unique_ranges.tsv"

# ----- Helpers (unchanged semantics) -----
def to_num(x: str):
    try:
        return float(x)
    except:
        return None

def to_int(x: str):
    try:
        return int(float(x))
    except:
        return None

def parse_tools(cell: str) -> List[str]:
    if not cell or cell == "NA":
        return []
    return [t.strip() for t in cell.split(",") if t.strip()]

# Structure per unique range:
# key: (start, end) -> info dict with flags and counts for classifying
# info = {
#   "partial_overlap": bool,               # participates in some partial overlap (<100)
#   "any_single_5tools_100": bool,         # seen as single-range with 5 tools and 100%
#   "single_tools_sets": list[set],        # observed tools in single-range rows
#   "tools_all_observed": set,             # all seen tools (from DF1 and/or DF2)
#   "promote_cat2": bool,                  # explicit promotion Cat3 -> Cat2 (rules below)
# }
def categorize(input_path: str, thresh: float, outdir: str = ".") -> None:
    os.makedirs(outdir, exist_ok=True)

    ranges: Dict[Tuple[int, int], Dict[str, object]] = {}

    with open(input_path, "r", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r, None)  # keep header if present
        for row in r:
            if len(row) < 12:
                continue

            # Columns per original script
            s1, e1, s2, e2 = row[0], row[1], row[2], row[3]
            ov1, ov2 = row[4], row[5]
            tools1, tools2 = row[10], row[11]

            # Parse numerics
            s1i, e1i = to_int(s1), to_int(e1)
            s2i, e2i = to_int(s2), to_int(e2)
            ov1f, ov2f = to_num(ov1), to_num(ov2)

            df1_valid = (s1 not in ("", "NA") and e1 not in ("", "NA") and s1i is not None and e1i is not None)
            df2_valid = (s2 not in ("", "NA") and e2 not in ("", "NA") and s2i is not None and e2i is not None)

            # Initialize and accumulate tools
            def ensure(key: Tuple[int, int]) -> None:
                if key not in ranges:
                    ranges[key] = {
                        "partial_overlap": False,
                        "any_single_5tools_100": False,
                        "single_tools_sets": [],
                        "tools_all_observed": set(),
                        "promote_cat2": False,
                    }

            if df1_valid:
                key1 = (s1i, e1i)
                ensure(key1)
                ranges[key1]["tools_all_observed"].update(parse_tools(tools1))

            if df2_valid:
                key2 = (s2i, e2i)
                ensure(key2)
                ranges[key2]["tools_all_observed"].update(parse_tools(tools2))

            # Row with two ranges (partial or total)
            if df1_valid and df2_valid:
                is_partial = (
                    (ov1f is not None and ov1f < 100.0) or
                    (ov2f is not None and ov2f < 100.0)
                )
                if is_partial:
                    ranges[(s1i, e1i)]["partial_overlap"] = True
                    ranges[(s2i, e2i)]["partial_overlap"] = True

                    # Promotion rules from Cat3 -> Cat2 according to threshold
                    n1 = len(parse_tools(tools1))
                    n2 = len(parse_tools(tools2))
                    both_have_overlap_values = (ov1f is not None and ov2f is not None)

                    # Rule A: 5 vs 5 and both overlaps <= THRESH -> PROMOTE both
                    if n1 == 5 and n2 == 5 and both_have_overlap_values:
                        if ov1f <= thresh and ov2f <= thresh:
                            ranges[(s1i, e1i)]["promote_cat2"] = True
                            ranges[(s2i, e2i)]["promote_cat2"] = True

                    # Rule B: 5 vs !=5 and both overlaps <= THRESH -> PROMOTE only the 5
                    elif (n1 == 5) ^ (n2 == 5) and both_have_overlap_values:
                        if ov1f <= thresh and ov2f <= thresh:
                            if n1 == 5:
                                ranges[(s1i, e1i)]["promote_cat2"] = True
                            else:
                                ranges[(s2i, e2i)]["promote_cat2"] = True

            # Rules with only one range (single-range)
            if df1_valid and not df2_valid:
                key1 = (s1i, e1i)
                t1 = set(parse_tools(tools1))
                ranges[key1]["single_tools_sets"].append(t1)
                if ov1f is not None and abs(ov1f - 100.0) < 1e-9 and len(t1) == 5:
                    ranges[key1]["any_single_5tools_100"] = True

            if df2_valid and not df1_valid:
                key2 = (s2i, e2i)
                t2 = set(parse_tools(tools2))
                ranges[key2]["single_tools_sets"].append(t2)
                if ov2f is not None and abs(ov2f - 100.0) < 1e-9 and len(t2) == 5:
                    ranges[key2]["any_single_5tools_100"] = True

    # ---------- Classification (priority preserved) ----------
    # 1) Cat2 if promote_cat2 == True  (explicit promotion)
    # 2) Cat3 if partial_overlap == True
    # 3) Cat2 if any_single_5tools_100 == True
    # 4) Cat1 if in all their single-range appearances had exactly 1 tool
    # 5) Cat4 the rest
    cat1: List[Tuple[int, int, int, str]] = []
    cat2: List[Tuple[int, int, int, str]] = []
    cat3: List[Tuple[int, int, int, str]] = []
    cat4: List[Tuple[int, int, int, str]] = []

    def write_tsv(path: str, rows: List[Tuple[int, int, int, str]]) -> None:
        with open(path, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            w.writerow(["Start", "End", "n_tools_observed", "tools_observed"])
            for (s, e, n, tools_str) in rows:
                w.writerow([s, e, n, tools_str])

    for (start, end), info in ranges.items():
        tools_all = sorted(info["tools_all_observed"])
        n_all = len(tools_all)
        tools_str = ",".join(tools_all)

        if info["promote_cat2"]:
            cat2.append((start, end, n_all, tools_str))
            continue

        if info["partial_overlap"]:
            cat3.append((start, end, n_all, tools_str))
            continue

        if info["any_single_5tools_100"]:
            cat2.append((start, end, n_all, tools_str))
            continue

        single_sets = info["single_tools_sets"]
        if single_sets and all(len(s) == 1 for s in single_sets):
            cat1.append((start, end, n_all, tools_str))
            continue

        cat4.append((start, end, n_all, tools_str))

    # ----- Write outputs (same filenames & columns) -----
    write_tsv(os.path.join(outdir, CAT1_FILE), cat1)
    write_tsv(os.path.join(outdir, CAT2_FILE), cat2)
    write_tsv(os.path.join(outdir, CAT3_FILE), cat3)
    write_tsv(os.path.join(outdir, CAT4_FILE), cat4)

    with open(os.path.join(outdir, ALL_UNIQUE_FILE), "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["Start", "End"])
        for (s, e) in sorted(ranges.keys()):
            w.writerow([s, e])

    # ----- Summary (same text) -----
    total = len(ranges)
    c1, c2n, c3, c4n = len(cat1), len(cat2), len(cat3), len(cat4)

    with open(os.path.join(outdir, SUMMARY_FILE), "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["category", "count"])
        w.writerow(["cat1_single_tool", c1])
        w.writerow(["cat2_full_consensus_5", c2n])
        w.writerow(["cat3_partial_overlaps", c3])
        w.writerow(["cat4_other_full_or_uncategorized", c4n])
        w.writerow(["TOTAL_unique_ranges", total])
        w.writerow(["CHECK_sum_categories", c1 + c2n + c3 + c4n])

    print("===== SUMMARY (by unique ranges) =====")
    print(f"Cat1 single-tool                  : {c1}")
    print(f"Cat2 full-consensus (5 tools)     : {c2n}")
    print(f"Cat3 partial overlaps (<100)      : {c3}")
    print(f"Cat4 other/full/uncategorized     : {c4n}")
    print("---------------------------------------")
    print(f"TOTAL unique ranges               : {total}")
    print(f"CHECK: sum categories             : {c1 + c2n + c3 + c4n}")
    print(f"(Promotions with THRESH={thresh}%) applied: 5vs5 -> both; 5vs≠5 -> only 5-tool range.")

def main():
    ap = argparse.ArgumentParser(description="Categorize unique ranges (original logic).")
    ap.add_argument("--input", default=DEFAULT_INPUT, help="Path to final_overlap_results_combined.tsv")
    ap.add_argument("--thresh", type=float, default=DEFAULT_THRESH, help="Promotion threshold in % (default 10.0)")
    ap.add_argument("--outdir", default=".", help="Output directory (default: current dir)")
    args = ap.parse_args()
    categorize(args.input, args.thresh, args.outdir)

if __name__ == "__main__":
    main()

