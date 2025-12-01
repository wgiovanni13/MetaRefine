#!/usr/bin/env python3
"""
Detect potential frameshifts by comparing adjacent genes' BLAST hits.

Strategy:
1. Extract genomic coordinates from query IDs
2. Identify adjacent genes (consecutive on genome)
3. Compare top 10 BLAST hits between adjacent genes
4. Find cases where adjacent genes hit the same reference protein
5. Check if they align to different, non-overlapping regions
6. Report potential frameshifts
"""

import argparse
import re
import sys
from collections import defaultdict

def parse_coordinates(query_id):
    """
    Extract genomic coordinates from query ID.
    
    Examples:
    - 176-545(-)_1 -> start=176, end=545, strand='-'
    - 687-1620(+)_1 -> start=687, end=1620, strand='+'
    
    Returns: (start, end, strand) or (None, None, None)
    """
    # Pattern: digits-digits(strand)
    match = re.match(r'(\d+)-(\d+)$$([+-])$$', query_id)
    
    if match:
        start = int(match.group(1))
        end = int(match.group(2))
        strand = match.group(3)
        return start, end, strand
    
    return None, None, None

def parse_blast_outfmt6(blast_file):
    """
    Parse BLAST output format 6 (tabular).
    
    Columns:
    0: qseqid, 1: sseqid, 2: pident, 3: length, 4: mismatch, 5: gapopen,
    6: qstart, 7: qend, 8: sstart, 9: send, 10: evalue, 11: bitscore,
    12: qlen, 13: slen, 14: stitle
    
    Returns: dict[query_id] = list of hit dicts
    """
    hits_by_query = defaultdict(list)
    
    with open(blast_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 15:
                continue
            
            query_id = fields[0]
            
            hit = {
                'qseqid': fields[0],
                'sseqid': fields[1],
                'pident': float(fields[2]),
                'length': int(fields[3]),
                'qstart': int(fields[6]),
                'qend': int(fields[7]),
                'sstart': int(fields[8]),
                'send': int(fields[9]),
                'evalue': float(fields[10]),
                'bitscore': float(fields[11]),
                'qlen': int(fields[12]),
                'slen': int(fields[13]),
                'stitle': fields[14] if len(fields) > 14 else ''
            }
            
            # Calculate coverage
            hit['qcov'] = (abs(hit['qend'] - hit['qstart']) + 1) / hit['qlen'] * 100
            hit['scov'] = (abs(hit['send'] - hit['sstart']) + 1) / hit['slen'] * 100
            
            hits_by_query[query_id].append(hit)
    
    return hits_by_query

def find_adjacent_genes(hits_by_query, max_gap=500):
    """
    Find pairs of adjacent genes on the genome.
    
    Adjacent = genes on the same strand with small gap between them.
    
    Returns: list of (gene1_id, gene2_id, gap_size) tuples
    """
    # Parse coordinates for all queries
    gene_positions = []
    
    for query_id in hits_by_query.keys():
        start, end, strand = parse_coordinates(query_id)
        
        if start is not None:
            gene_positions.append({
                'id': query_id,
                'start': start,
                'end': end,
                'strand': strand
            })
    
    # Sort by strand, then by start position
    gene_positions.sort(key=lambda x: (x['strand'], x['start']))
    
    # Find adjacent pairs
    adjacent_pairs = []
    
    for i in range(len(gene_positions) - 1):
        gene1 = gene_positions[i]
        gene2 = gene_positions[i + 1]
        
        # Must be on same strand
        if gene1['strand'] != gene2['strand']:
            continue
        
        # Calculate gap
        if gene1['strand'] == '+':
            gap = gene2['start'] - gene1['end']
        else:
            # For minus strand, "adjacent" means gene2 starts where gene1 ends
            gap = gene2['start'] - gene1['end']
        
        # Filter by max gap
        if gap >= 0 and gap <= max_gap:
            adjacent_pairs.append((gene1['id'], gene2['id'], gap))
    
    return adjacent_pairs

def detect_frameshifts(hits_by_query, adjacent_pairs, 
                       min_identity=30.0,
                       min_coverage=50.0,
                       max_overlap=50):
    """
    Detect potential frameshifts in adjacent gene pairs.
    
    Criteria:
    1. Both genes hit the same reference protein
    2. Hits have sufficient identity and coverage
    3. Alignments cover different regions of the reference (minimal overlap)
    4. Combined coverage suggests they're fragments of the same gene
    """
    frameshifts = []
    
    for gene1_id, gene2_id, gap in adjacent_pairs:
        # Get top hits for each gene
        hits1 = hits_by_query.get(gene1_id, [])
        hits2 = hits_by_query.get(gene2_id, [])
        
        if not hits1 or not hits2:
            continue
        
        # Compare all hits
        for hit1 in hits1[:10]:  # Top 10 hits
            for hit2 in hits2[:10]:
                # Check if same reference protein
                if hit1['sseqid'] != hit2['sseqid']:
                    continue
                
                # Check identity thresholds
                if hit1['pident'] < min_identity or hit2['pident'] < min_identity:
                    continue
                
                # Check coverage thresholds
                if hit1['qcov'] < min_coverage or hit2['qcov'] < min_coverage:
                    continue
                
                # Get alignment positions on reference
                ref_start1 = min(hit1['sstart'], hit1['send'])
                ref_end1 = max(hit1['sstart'], hit1['send'])
                ref_start2 = min(hit2['sstart'], hit2['send'])
                ref_end2 = max(hit2['sstart'], hit2['send'])
                
                # Calculate overlap on reference
                overlap_start = max(ref_start1, ref_start2)
                overlap_end = min(ref_end1, ref_end2)
                overlap_len = max(0, overlap_end - overlap_start + 1)
                
                # Check if minimal overlap (should be different regions)
                if overlap_len > max_overlap:
                    continue
                
                # Calculate combined coverage
                combined_start = min(ref_start1, ref_start2)
                combined_end = max(ref_end1, ref_end2)
                combined_coverage = (combined_end - combined_start + 1) / hit1['slen'] * 100
                
                # Store frameshift candidate
                frameshifts.append({
                    'gene1_id': gene1_id,
                    'gene2_id': gene2_id,
                    'gap_nt': gap,
                    'ref_id': hit1['sseqid'],
                    'ref_length': hit1['slen'],
                    'ref_description': hit1['stitle'],
                    'gene1_ref_start': ref_start1,
                    'gene1_ref_end': ref_end1,
                    'gene1_identity': hit1['pident'],
                    'gene1_qcov': hit1['qcov'],
                    'gene1_scov': hit1['scov'],
                    'gene2_ref_start': ref_start2,
                    'gene2_ref_end': ref_end2,
                    'gene2_identity': hit2['pident'],
                    'gene2_qcov': hit2['qcov'],
                    'gene2_scov': hit2['scov'],
                    'overlap_aa': overlap_len,
                    'combined_coverage': combined_coverage,
                    'gene1_evalue': hit1['evalue'],
                    'gene2_evalue': hit2['evalue'],
                    'gene1_bitscore': hit1['bitscore'],
                    'gene2_bitscore': hit2['bitscore']
                })
    
    return frameshifts

def main():
    parser = argparse.ArgumentParser(
        description='Detect frameshifts from BLAST results against NR',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python3 detect_frameshifts_from_nr.py \\
      --blast proteins_vs_nr.blastp \\
      --output frameshifts_nr.tsv \\
      --max-gap 500 \\
      --min-identity 30.0 \\
      --min-coverage 50.0
        """
    )
    
    parser.add_argument('--blast', required=True,
                       help='BLAST output file (format 6)')
    parser.add_argument('--output', default='frameshifts_nr.tsv',
                       help='Output TSV file')
    parser.add_argument('--max-gap', type=int, default=500,
                       help='Max gap (nt) between adjacent genes (default: 500)')
    parser.add_argument('--min-identity', type=float, default=30.0,
                       help='Min percent identity (default: 30.0)')
    parser.add_argument('--min-coverage', type=float, default=50.0,
                       help='Min query coverage (default: 50.0)')
    parser.add_argument('--max-overlap', type=int, default=50,
                       help='Max overlap (aa) on reference (default: 50)')
    
    args = parser.parse_args()
    
    print("=== FRAMESHIFT DETECTION FROM NR BLAST ===")
    print(f"BLAST file: {args.blast}")
    print(f"Max gap between genes: {args.max_gap} nt")
    print(f"Min identity: {args.min_identity}%")
    print(f"Min coverage: {args.min_coverage}%")
    print(f"Max overlap on reference: {args.max_overlap} aa")
    print()
    
    # Parse BLAST results
    print("Parsing BLAST results...")
    hits_by_query = parse_blast_outfmt6(args.blast)
    print(f"Queries with hits: {len(hits_by_query):,}")
    
    total_hits = sum(len(hits) for hits in hits_by_query.values())
    print(f"Total hits: {total_hits:,}")
    print()
    
    # Find adjacent genes
    print("Finding adjacent genes...")
    adjacent_pairs = find_adjacent_genes(hits_by_query, args.max_gap)
    print(f"Adjacent gene pairs: {len(adjacent_pairs):,}")
    print()
    
    # Detect frameshifts
    print("Detecting frameshifts...")
    frameshifts = detect_frameshifts(
        hits_by_query,
        adjacent_pairs,
        args.min_identity,
        args.min_coverage,
        args.max_overlap
    )
    
    print(f"Potential frameshifts detected: {len(frameshifts):,}")
    print()
    
    # Write output
    if frameshifts:
        with open(args.output, 'w') as f:
            # Header
            f.write('\t'.join([
                'gene1_id', 'gene2_id', 'gap_nt',
                'ref_id', 'ref_length', 'ref_description',
                'gene1_ref_start', 'gene1_ref_end', 'gene1_identity', 'gene1_qcov', 'gene1_scov',
                'gene2_ref_start', 'gene2_ref_end', 'gene2_identity', 'gene2_qcov', 'gene2_scov',
                'overlap_aa', 'combined_coverage',
                'gene1_evalue', 'gene2_evalue', 'gene1_bitscore', 'gene2_bitscore'
            ]) + '\n')
            
            # Data
            for fs in frameshifts:
                f.write('\t'.join([
                    fs['gene1_id'], fs['gene2_id'], str(fs['gap_nt']),
                    fs['ref_id'], str(fs['ref_length']), fs['ref_description'],
                    str(fs['gene1_ref_start']), str(fs['gene1_ref_end']),
                    f"{fs['gene1_identity']:.2f}", f"{fs['gene1_qcov']:.2f}", f"{fs['gene1_scov']:.2f}",
                    str(fs['gene2_ref_start']), str(fs['gene2_ref_end']),
                    f"{fs['gene2_identity']:.2f}", f"{fs['gene2_qcov']:.2f}", f"{fs['gene2_scov']:.2f}",
                    str(fs['overlap_aa']), f"{fs['combined_coverage']:.2f}",
                    f"{fs['gene1_evalue']:.2e}", f"{fs['gene2_evalue']:.2e}",
                    f"{fs['gene1_bitscore']:.1f}", f"{fs['gene2_bitscore']:.1f}"
                ]) + '\n')
        
        print(f"Output written to: {args.output}")
        
        # Show examples
        print("\n=== TOP 10 FRAMESHIFT CANDIDATES ===")
        print(f"{'Gene1':<20} {'Gene2':<20} {'Gap':<6} {'RefID':<15} {'Combined_Cov':<12} {'Description':<50}")
        print("-" * 130)
        
        # Sort by combined coverage (descending)
        sorted_fs = sorted(frameshifts, key=lambda x: x['combined_coverage'], reverse=True)
        
        for fs in sorted_fs[:10]:
            desc = fs['ref_description'][:47] + '...' if len(fs['ref_description']) > 50 else fs['ref_description']
            print(f"{fs['gene1_id']:<20} {fs['gene2_id']:<20} {fs['gap_nt']:<6} "
                  f"{fs['ref_id']:<15} {fs['combined_coverage']:>6.1f}%     {desc}")
    
    else:
        print("No frameshifts detected!")
        # Create empty output file
        with open(args.output, 'w') as f:
            f.write('# No frameshifts detected\n')

if __name__ == '__main__':
    main()

