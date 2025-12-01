#!/usr/bin/env python3
"""
Select best query per group based on HHpred results.
Expected format: start-end_(strand)_frame.hhr
Example: 2730515-2731471_(-)_4.hhr
"""
import os, re, math, sys, argparse
from collections import defaultdict

# Regex para nombres: start-end_(+/-)_frame.hhr
NAME_RE = re.compile(r'^(\d+)-(\d+)_\(([+-])\)_(\d+)\.hhr$')

# Detecta el encabezado de la tabla de HHpred
HEADER_RE = re.compile(r'^\s*No\b.*\bHit\b.*\bProb\b.*\bE-value\b.*\bP-value\b.*\bScore\b', re.I)

# Fila de resultado del top table
ROW_RE = re.compile(
    r'^\s*\d+\s+(\S+)\s+.+?\s+([0-9.]+)\s+([0-9.eE+\-]+)\s+([0-9.eE+\-]+)\s+([0-9.]+)\s+([0-9.]+)\s+(\d+)\s+(\d+)-(\d+)',
    re.ASCII
)

def parse_args():
    parser = argparse.ArgumentParser(
        description='Select best query per group based on HHpred results'
    )
    parser.add_argument('--groups', required=True, help='Path to groups.tsv file')
    parser.add_argument('--hhrdir', required=True, help='Directory containing .hhr files')
    parser.add_argument('--outdir', default='.', help='Output directory (default: current directory)')
    parser.add_argument('--max-hits', type=int, default=10, help='Maximum number of hits to consider (default: 10)')
    return parser.parse_args()

def parse_groups(groups_path):
    """Lee groups.tsv y devuelve dict[group_id] -> list[basename.hhr]"""
    groups = {}
    with open(groups_path, "r") as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            gid = parts[0]
            files_field = parts[6]
            files = [x.strip() for x in files_field.split(";") if x.strip()]
            groups[gid] = files
    return groups

def parse_top_hits(hhr_path, max_hits=10):
    """Devuelve lista de dicts con hits de HHpred"""
    hits = []
    in_table = False
    with open(hhr_path, 'r', errors="ignore") as f:
        for line in f:
            if not in_table:
                if HEADER_RE.search(line):
                    in_table = True
                continue
            if len(hits) >= max_hits:
                break
            m = ROW_RE.match(line)
            if not m:
                continue
            hit_id = m.group(1)
            try:
                prob = float(m.group(2))
                evalue = float(m.group(3))
                pvalue = float(m.group(4))
                score = float(m.group(5))
                cols = int(m.group(7))
                qstart = int(m.group(8))
                qend = int(m.group(9))
            except Exception:
                continue
            hits.append({
                "id": hit_id,
                "prob": prob,
                "evalue": evalue,
                "pvalue": pvalue,
                "score": score,
                "cols": cols,
                "qstart": qstart,
                "qend": qend
            })
    return hits

def choose_best_in_group(basenames, hhrdir, max_hits=10):
    """Selecciona el mejor query de un grupo"""
    enriched = []
    for base in basenames:
        m = NAME_RE.match(base)
        if not m:
            continue
        start = int(m.group(1))
        end = int(m.group(2))
        strand = m.group(3)
        frame = int(m.group(4))
        
        path = os.path.join(hhrdir, base)
        if not os.path.exists(path):
            continue
        
        hits = parse_top_hits(path, max_hits)
        if not hits:
            continue
        
        enriched.append({
            "file": base,
            "start": start,
            "end": end,
            "strand": strand,
            "frame": frame,
            "top1": hits[0],
            "top_hits": hits
        })
    
    if not enriched:
        return None

    best_prob = max(x["top1"]["prob"] for x in enriched)
    best_score = max(x["top1"]["score"] for x in enriched)
    best_e = min(x["top1"]["evalue"] for x in enriched)
    best_p = min(x["top1"]["pvalue"] for x in enriched)

    for x in enriched:
        t = x["top1"]
        wins = 0
        if abs(t["prob"] - best_prob) < 1e-9:
            wins += 1
        if abs(t["evalue"] - best_e) < 1e-300 or t["evalue"] == best_e:
            wins += 1
        if abs(t["pvalue"] - best_p) < 1e-300 or t["pvalue"] == best_p:
            wins += 1
        if abs(t["score"] - best_score) < 1e-9:
            wins += 1
        x["wins4"] = wins

    n = len(enriched)
    thresh = math.ceil((n - 1) / 2.0)
    
    for i, x in enumerate(enriched):
        t1_id = x["top1"]["id"]
        support = 0
        for j, y in enumerate(enriched):
            if i == j:
                continue
            if any(h["id"] == t1_id for h in y["top_hits"]):
                support += 1
        x["support_in_others"] = support
        x["support_thresh"] = thresh
        x["has_majority_support"] = (support >= thresh)

    def sort_key(x):
        t = x["top1"]
        return (-x["wins4"], -t["prob"], t["evalue"], t["pvalue"], -t["score"], -t["cols"])

    cand_majority = [x for x in enriched if x["has_majority_support"]]
    if cand_majority:
        best = sorted(cand_majority, key=sort_key)[0]
    else:
        cand_wins = [x for x in enriched if x["wins4"] >= 3]
        if cand_wins:
            best = sorted(cand_wins, key=sort_key)[0]
        else:
            best = sorted(enriched, key=sort_key)[0]

    t = best["top1"]
    high_conf = 1 if (best["has_majority_support"] and best["wins4"] >= 3) else 0

    return {
        "best_file": best["file"],
        "strand": best["strand"],
        "start": best["start"],
        "end": best["end"],
        "frame": best["frame"],
        "top1_id": t["id"],
        "prob": t["prob"],
        "evalue": t["evalue"],
        "pvalue": t["pvalue"],
        "score": t["score"],
        "cols": t["cols"],
        "wins4": best["wins4"],
        "support_in_others": best["support_in_others"],
        "support_thresh": best["support_thresh"],
        "group_size": n,
        "high_confidence": high_conf
    }

def main():
    args = parse_args()
    
    if not os.path.isfile(args.groups):
        print(f"ERROR: Groups file not found: {args.groups}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.isdir(args.hhrdir):
        print(f"ERROR: HHR directory not found: {args.hhrdir}", file=sys.stderr)
        sys.exit(1)
    
    os.makedirs(args.outdir, exist_ok=True)
    
    groups = parse_groups(args.groups)
    print(f"Loaded {len(groups)} groups from {args.groups}")
    
    out = []
    skipped_groups = 0
    
    for gid, files in groups.items():
        info = choose_best_in_group(files, args.hhrdir, args.max_hits)
        if info is None:
            skipped_groups += 1
            continue
        info["group_id"] = gid
        out.append(info)

    out_path = os.path.join(args.outdir, "best_query_per_group.tsv")
    hi = sum(1 for r in out if r["high_confidence"] == 1)
    
    with open(out_path, "w") as f:
        header = ["group_id", "group_size", "strand", "best_query", "start", "end", "frame",
                  "top1_id", "prob", "evalue", "pvalue", "score", "cols",
                  "wins4", "support_in_others", "support_thresh", "high_confidence"]
        f.write("\t".join(header) + "\n")
        
        def gid_key(x):
            try:
                return int(x["group_id"])
            except:
                return x["group_id"]
        
        for r in sorted(out, key=gid_key):
            row = [
                str(r["group_id"]), str(r["group_size"]), r["strand"],
                r["best_file"], str(r["start"]), str(r["end"]), str(r["frame"]),
                r["top1_id"],
                f'{r["prob"]:.1f}', f'{r["evalue"]:.3e}', f'{r["pvalue"]:.3e}',
                f'{r["score"]:.1f}', str(r["cols"]),
                str(r["wins4"]), str(r["support_in_others"]),
                str(r["support_thresh"]), str(r["high_confidence"])
            ]
            f.write("\t".join(row) + "\n")

    print(f"Done: {out_path}")
    print(f"Groups evaluated     : {len(out)}")
    print(f"Groups without data  : {skipped_groups}")
    print(f"High confidence      : {hi} / {len(out)}")

if __name__ == "__main__":
    main()

