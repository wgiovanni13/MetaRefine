#!/usr/bin/env python3
"""
Frameshift detector from DIAMOND results
Identifies adjacent genes that align to the same reference protein
"""

import sys
import gzip
from collections import defaultdict

def open_file(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    return open(filename, 'r')

def parse_diamond_file(diamond_file):
    """Parse DIAMOND tabular output and store gene hits"""
    gene_hits = defaultdict(list)
    with open_file(diamond_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 15:  # Need 15 columns (including description)
                continue
            
            gene_hits[fields[0]].append({
                'target_id': fields[1],
                'target_desc': fields[14] if len(fields) > 14 else '',
                'pident': float(fields[2]),
                't_start': int(fields[8]),
                't_end': int(fields[9]),
                'evalue': float(fields[10]),
                't_len': int(fields[13])
            })
    
    print(f"[INFO] Loaded {len(gene_hits)} genes", file=sys.stderr)
    return gene_hits

def extract_gene_coords(gene_id):
    """Extract genomic coordinates from gene ID (e.g., 678032-678770(-)_1)"""
    try:
        parts = gene_id.split('-', 1)
        if len(parts) != 2:
            return None, None, None, None
        
        start = int(parts[0])
        rest = parts[1]
        paren_pos = rest.find('(')
        if paren_pos == -1:
            return None, None, None, None
        
        end = int(rest[:paren_pos])
        strand = rest[paren_pos + 1]
        underscore_pos = rest.find('_')
        if underscore_pos == -1:
            return None, None, None, None
        
        frame = int(rest[underscore_pos + 1:])
        return start, end, strand, frame
    except (ValueError, IndexError):
        return None, None, None, None

def is_characterized(description):
    """Check if protein is characterized (not hypothetical/uncharacterized)"""
    desc_lower = description.lower()
    uncharacterized_terms = [
        'uncharacterized',
        'hypothetical',
        'unknown function',
        'putative uncharacterized',
        'duf',
        'predicted protein'
    ]
    return not any(term in desc_lower for term in uncharacterized_terms)

def calculate_overlap(hit1, hit2):
    """Calculate overlap fraction between two alignments on reference"""
    overlap_start = max(hit1['t_start'], hit2['t_start'])
    overlap_end = min(hit1['t_end'], hit2['t_end'])
    overlap = max(0, overlap_end - overlap_start + 1)
    len1 = hit1['t_end'] - hit1['t_start'] + 1
    len2 = hit2['t_end'] - hit2['t_start'] + 1
    return overlap/len1 if len1 > 0 else 0, overlap/len2 if len2 > 0 else 0, overlap

def score_frameshift(fs):
    """Calculate quality score for frameshift candidate"""
    # Prioritize: characterized > avg_identity > total_cov
    char_bonus = 1000 if is_characterized(fs['target_desc']) else 0
    return char_bonus + fs['avg_identity'] + fs['total_cov']

def find_frameshifts(gene_hits, top_n=20, min_identity=30, max_evalue=1e-3, max_overlap=0.3):
    """Find potential frameshifts in adjacent genes"""
    
    # Sort genes by genomic position
    sorted_genes = sorted(gene_hits.keys(), key=lambda g: extract_gene_coords(g)[0] or float('inf'))
    valid_genes = [g for g in sorted_genes if extract_gene_coords(g)[0] is not None]
    
    print(f"[INFO] Analyzing {len(valid_genes)} valid genes", file=sys.stderr)
    
    # Store best frameshift per gene pair
    pair_frameshifts = {}
    
    for i in range(len(valid_genes) - 1):
        gene1 = valid_genes[i]
        gene2 = valid_genes[i+1]
        
        s1, e1, strand1, _ = extract_gene_coords(gene1)
        s2, e2, strand2, _ = extract_gene_coords(gene2)
        
        # Filter and limit hits
        hits1 = [h for h in gene_hits[gene1] if h['pident'] >= min_identity and h['evalue'] <= max_evalue][:top_n]
        hits2 = [h for h in gene_hits[gene2] if h['pident'] >= min_identity and h['evalue'] <= max_evalue][:top_n]
        
        if not hits1 or not hits2:
            continue
        
        # Find shared targets
        candidates = []
        for hit1 in hits1:
            for hit2 in hits2:
                if hit1['target_id'] != hit2['target_id']:
                    continue
                
                ovlp1, ovlp2, ovlp_len = calculate_overlap(hit1, hit2)
                
                # Filter by overlap
                if ovlp1 > max_overlap and ovlp2 > max_overlap:
                    continue
                
                # Calculate metrics
                ref_len = hit1['t_len']
                cov1 = (hit1['t_end'] - hit1['t_start'] + 1) / ref_len * 100
                cov2 = (hit2['t_end'] - hit2['t_start'] + 1) / ref_len * 100
                ref_start = min(hit1['t_start'], hit2['t_start'])
                ref_end = max(hit1['t_end'], hit2['t_end'])
                total_cov = (ref_end - ref_start + 1) / ref_len * 100
                
                gap = (hit2['t_start'] - hit1['t_end'] - 1) if hit1['t_start'] < hit2['t_start'] else (hit1['t_start'] - hit2['t_end'] - 1)
                ref_order = "1->2" if hit1['t_start'] < hit2['t_start'] else "2->1"
                avg_identity = (hit1['pident'] + hit2['pident']) / 2
                
                if total_cov >= 70 and avg_identity >= 50:
                    confidence = "HIGH"
                elif total_cov >= 50 and avg_identity >= 40:
                    confidence = "MEDIUM"
                else:
                    confidence = "LOW"
                
                candidates.append({
                    'gene1': gene1, 'gene2': gene2,
                    'gene1_coords': f"{s1}-{e1}", 'gene2_coords': f"{s2}-{e2}",
                    'gene1_strand': strand1, 'gene2_strand': strand2,
                    'distance_bp': s2 - e1 - 1,
                    'ref_id': hit1['target_id'],
                    'target_desc': hit1['target_desc'],
                    'ref_length': ref_len,
                    'gene1_identity': hit1['pident'], 'gene2_identity': hit2['pident'], 'avg_identity': avg_identity,
                    'gene1_evalue': hit1['evalue'], 'gene2_evalue': hit2['evalue'],
                    'gene1_ref_start': hit1['t_start'], 'gene1_ref_end': hit1['t_end'],
                    'gene2_ref_start': hit2['t_start'], 'gene2_ref_end': hit2['t_end'],
                    'gene1_cov': cov1, 'gene2_cov': cov2, 'total_cov': total_cov,
                    'gap_in_ref': gap, 'overlap_aa': ovlp_len,
                    'ref_order': ref_order, 'confidence': confidence
                })
        
        if candidates:
            # Sort by score (characterized proteins get +1000 bonus)
            candidates.sort(key=score_frameshift, reverse=True)
            
            # ONLY accept characterized proteins
            best_characterized = None
            for candidate in candidates:
                if is_characterized(candidate['target_desc']):
                    best_characterized = candidate
                    break
            
            # Only add to results if we found a characterized protein
            if best_characterized:
                pair_key = (gene1, gene2)
                pair_frameshifts[pair_key] = best_characterized
    
    frameshifts = list(pair_frameshifts.values())
    print(f"[INFO] Found {len(frameshifts)} frameshift candidates", file=sys.stderr)
    return frameshifts

def write_tsv(frameshifts, output_file):
    """Write results to TSV file"""
    with open(output_file, 'w') as f:
        f.write("gene1\tgene2\tgene1_coords\tgene2_coords\tgene1_strand\tgene2_strand\tdistance_bp\t"
               "ref_id\ttarget_description\tref_length\t"
               "gene1_identity\tgene2_identity\tavg_identity\t"
               "gene1_evalue\tgene2_evalue\t"
               "gene1_ref_start\tgene1_ref_end\tgene2_ref_start\tgene2_ref_end\t"
               "gene1_cov\tgene2_cov\ttotal_cov\t"
               "gap_in_ref\toverlap_aa\tref_order\tconfidence\n")
        
        for fs in frameshifts:
            f.write(f"{fs['gene1']}\t{fs['gene2']}\t{fs['gene1_coords']}\t{fs['gene2_coords']}\t"
                   f"{fs['gene1_strand']}\t{fs['gene2_strand']}\t{fs['distance_bp']}\t"
                   f"{fs['ref_id']}\t{fs['target_desc']}\t{fs['ref_length']}\t"
                   f"{fs['gene1_identity']:.1f}\t{fs['gene2_identity']:.1f}\t{fs['avg_identity']:.1f}\t"
                   f"{fs['gene1_evalue']:.2e}\t{fs['gene2_evalue']:.2e}\t"
                   f"{fs['gene1_ref_start']}\t{fs['gene1_ref_end']}\t{fs['gene2_ref_start']}\t{fs['gene2_ref_end']}\t"
                   f"{fs['gene1_cov']:.1f}\t{fs['gene2_cov']:.1f}\t{fs['total_cov']:.1f}\t"
                   f"{fs['gap_in_ref']}\t{fs['overlap_aa']}\t{fs['ref_order']}\t{fs['confidence']}\n")

def main():
    if len(sys.argv) < 3:
        print("Usage: detect_frameshifts.py <input.diamond.gz> <output.tsv>", file=sys.stderr)
        sys.exit(1)
    
    gene_hits = parse_diamond_file(sys.argv[1])
    
    if len(gene_hits) == 0:
        # Write empty file with header
        with open(sys.argv[2], 'w') as f:
            f.write("gene1\tgene2\tgene1_coords\tgene2_coords\tgene1_strand\tgene2_strand\tdistance_bp\t"
                   "ref_id\ttarget_description\tref_length\t"
                   "gene1_identity\tgene2_identity\tavg_identity\t"
                   "gene1_evalue\tgene2_evalue\t"
                   "gene1_ref_start\tgene1_ref_end\tgene2_ref_start\tgene2_ref_end\t"
                   "gene1_cov\tgene2_cov\ttotal_cov\t"
                   "gap_in_ref\toverlap_aa\tref_order\tconfidence\n")
        return
    
    frameshifts = find_frameshifts(gene_hits)
    write_tsv(frameshifts, sys.argv[2])

if __name__ == '__main__':
    main()

