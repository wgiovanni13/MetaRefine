#!/usr/bin/env python3

"""
Detect potential frameshifts from BLAST all-vs-all results.
Logic: Find cases where one sequence (query) aligns to two fragments (subjects)
that are close together and together cover most of the query sequence.
The query should be LONGER than individual subjects (query is the complete gene).
"""

import sys
import argparse
from collections import defaultdict
import re

def parse_range_id(range_id):
    """
    Parse range ID to extract coordinates.
    Handles multiple formats:
    - Format 1: start-end_strand_frame::contig:start-end(strand)
    - Format 2: contig:start-end(strand)
    - Format 3: start-end(strand)_frame
    Returns: (start, end, strand, contig)
    """
    # Try format: 792151-792423_+_1::contig_1_np1212_8_pilon:792150-792423(+)
    # Extract from the part after ::
    if '::' in range_id:
        parts = range_id.split('::')
        if len(parts) == 2:
            # Parse the second part: contig:start-end(strand)
            match = re.match(r'(.+):(\d+)-(\d+)\(([+-])\)', parts[1])
            if match:
                contig = match.group(1)
                start = int(match.group(2))
                end = int(match.group(3))
                strand = match.group(4)
                return start, end, strand, contig
    
    # Try format: contig:start-end(strand)
    match = re.match(r'(.+):(\d+)-(\d+)\(([+-])\)', range_id)
    if match:
        contig = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        strand = match.group(4)
        return start, end, strand, contig
    
    # Try format: start-end(strand)_frame
    match = re.match(r'(\d+)-(\d+)\(([+-])\)_(\d+)', range_id)
    if match:
        start = int(match.group(1))
        end = int(match.group(2))
        strand = match.group(3)
        return start, end, strand, "unknown"
    
    # Try format: start-end_strand_frame
    match = re.match(r'(\d+)-(\d+)_([+-])_(\d+)', range_id)
    if match:
        start = int(match.group(1))
        end = int(match.group(2))
        strand = match.group(3)
        return start, end, strand, "unknown"
    
    return None, None, None, None

def load_blast_hits(blast_file):
    """
    Load BLAST hits from filtered file.
    Returns dict: query -> list of subject hits
    """
    hits_by_query = defaultdict(list)
    
    with open(blast_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 14:
                continue
            
            qseqid = fields[0]
            sseqid = fields[1]
            pident = float(fields[2])
            length = int(fields[3])
            qstart = int(fields[6])
            qend = int(fields[7])
            sstart = int(fields[8])
            send = int(fields[9])
            evalue = float(fields[10])
            bitscore = float(fields[11])
            qlen = int(fields[12])
            slen = int(fields[13])
            
            # Skip self-hits (should already be filtered, but double-check)
            if qseqid == sseqid:
                continue
            
            hit = {
                'sseqid': sseqid,
                'pident': pident,
                'length': length,
                'qstart': qstart,
                'qend': qend,
                'sstart': sstart,
                'send': send,
                'evalue': evalue,
                'bitscore': bitscore,
                'qlen': qlen,
                'slen': slen
            }
            
            hits_by_query[qseqid].append(hit)
    
    return hits_by_query

def calculate_genomic_gap(s1_start, s1_end, s2_start, s2_end):
    """
    Calculate genomic gap between two ranges.
    Returns gap size (negative if overlapping).
    """
    if s2_start > s1_end:
        return s2_start - s1_end
    elif s1_start > s2_end:
        return s1_start - s2_end
    else:
        # Overlapping - return negative value
        return min(s1_end, s2_end) - max(s1_start, s2_start)

def calculate_query_coverage(qstart1, qend1, qstart2, qend2):
    """
    Calculate combined coverage of query by two alignments.
    Returns: (coverage_length, overlap, gap_in_query)
    """
    # Ensure start < end
    if qstart1 > qend1:
        qstart1, qend1 = qend1, qstart1
    if qstart2 > qend2:
        qstart2, qend2 = qend2, qstart2
    
    # Calculate overlap in query
    overlap = 0
    if qstart2 <= qend1 and qstart1 <= qend2:
        overlap = min(qend1, qend2) - max(qstart1, qstart2)
    
    # Calculate gap in query
    gap_in_query = 0
    if qend1 < qstart2:
        gap_in_query = qstart2 - qend1
    elif qend2 < qstart1:
        gap_in_query = qstart1 - qend2
    
    # Total coverage (without double-counting overlap)
    coverage = (qend1 - qstart1) + (qend2 - qstart2) - overlap
    
    return coverage, overlap, gap_in_query

def detect_frameshifts(hits_by_query, min_identity, max_gap, min_query_coverage, max_subject_ratio):
    """
    Detect frameshift candidates.
    
    Criteria:
    1. One query aligns to two subjects
    2. Query should be longer than individual subjects (query = complete gene)
    3. Subjects are on same contig and strand
    4. Subjects are separated by a small gap (or overlap slightly)
    5. Combined coverage of query is high (min_query_coverage %)
    6. Alignments cover different parts of the query
    7. Each subject should not be much longer than the query
    """
    candidates = []
    
    for query, hits in hits_by_query.items():
        # Need at least 2 hits to detect frameshift
        if len(hits) < 2:
            continue
        
        # Parse query coordinates
        q_start, q_end, q_strand, q_contig = parse_range_id(query)
        if q_start is None:
            continue
        
        q_len = q_end - q_start
        
        # Compare all pairs of subjects
        for i, hit1 in enumerate(hits):
            # Filter by identity
            if hit1['pident'] < min_identity:
                continue
            
            # Parse subject 1 coordinates
            s1_start, s1_end, s1_strand, s1_contig = parse_range_id(hit1['sseqid'])
            if s1_start is None:
                continue
            
            s1_len = s1_end - s1_start
            
            # CRITICAL: Subject should be shorter than query
            # (query is the complete gene, subject is a fragment)
            if s1_len > q_len * max_subject_ratio:
                continue
            
            for hit2 in hits[i+1:]:
                # Filter by identity
                if hit2['pident'] < min_identity:
                    continue
                
                # Parse subject 2 coordinates
                s2_start, s2_end, s2_strand, s2_contig = parse_range_id(hit2['sseqid'])
                if s2_start is None:
                    continue
                
                s2_len = s2_end - s2_start
                
                # CRITICAL: Subject should be shorter than query
                if s2_len > q_len * max_subject_ratio:
                    continue
                
                # Must be on same contig
                if s1_contig != s2_contig:
                    continue
                
                # Must be on same strand
                if s1_strand != s2_strand:
                    continue
                
                # Calculate genomic gap between subjects
                genomic_gap = calculate_genomic_gap(s1_start, s1_end, s2_start, s2_end)
                
                # Check gap criteria (allow small overlap or gap up to max_gap)
                if genomic_gap > max_gap:
                    continue
                
                # Calculate query coverage
                query_coverage, query_overlap, query_gap = calculate_query_coverage(
                    hit1['qstart'], hit1['qend'],
                    hit2['qstart'], hit2['qend']
                )
                
                # Calculate coverage percentage
                coverage_pct = (query_coverage / hit1['qlen']) * 100
                
                # Check if combined coverage is sufficient
                if coverage_pct < min_query_coverage:
                    continue
                
                # Check that alignments cover different parts of query (not too much overlap)
                if query_overlap > 0.5 * min(hit1['length'], hit2['length']):
                    continue  # Too much overlap in query coverage
                
                # Found a frameshift candidate!
                candidates.append({
                    'query': query,
                    'q_contig': q_contig,
                    'q_start': q_start,
                    'q_end': q_end,
                    'q_len': q_len,
                    'q_strand': q_strand,
                    'subject1': hit1['sseqid'],
                    's1_contig': s1_contig,
                    's1_start': s1_start,
                    's1_end': s1_end,
                    's1_len': s1_len,
                    's1_identity': hit1['pident'],
                    's1_qstart': hit1['qstart'],
                    's1_qend': hit1['qend'],
                    'subject2': hit2['sseqid'],
                    's2_contig': s2_contig,
                    's2_start': s2_start,
                    's2_end': s2_end,
                    's2_len': s2_len,
                    's2_identity': hit2['pident'],
                    's2_qstart': hit2['qstart'],
                    's2_qend': hit2['qend'],
                    'genomic_gap': genomic_gap,
                    'query_coverage': query_coverage,
                    'query_coverage_pct': coverage_pct,
                    'query_overlap': query_overlap,
                    'query_gap': query_gap,
                    'combined_len': s1_len + s2_len,
                    'strand': s1_strand
                })
    
    return candidates

def parse_args():
    parser = argparse.ArgumentParser(
        description='Detect frameshift candidates from BLAST results'
    )
    parser.add_argument('--blast', required=True,
                       help='Filtered BLAST output file')
    parser.add_argument('--output', default='frameshift_candidates.txt',
                       help='Output text file')
    parser.add_argument('--output-tsv', default='frameshift_candidates.tsv',
                       help='Output TSV file')
    parser.add_argument('--min-identity', type=float, default=70.0,
                       help='Minimum percent identity (default: 70)')
    parser.add_argument('--max-gap', type=int, default=200,
                       help='Maximum genomic gap between fragments (default: 200 bp)')
    parser.add_argument('--min-query-coverage', type=float, default=80.0,
                       help='Minimum query coverage by both fragments (default: 80%%)')
    parser.add_argument('--max-subject-ratio', type=float, default=1.2,
                       help='Maximum ratio of subject length to query length (default: 1.2)')
    return parser.parse_args()

def main():
    args = parse_args()
    
    print(f"Loading BLAST results from: {args.blast}")
    hits_by_query = load_blast_hits(args.blast)
    print(f"  Loaded hits for {len(hits_by_query)} queries")
    
    print(f"\nDetecting frameshifts...")
    print(f"  Min identity: {args.min_identity}%")
    print(f"  Max genomic gap: {args.max_gap} bp")
    print(f"  Min query coverage: {args.min_query_coverage}%")
    print(f"  Max subject/query ratio: {args.max_subject_ratio}")
    
    candidates = detect_frameshifts(
        hits_by_query,
        args.min_identity,
        args.max_gap,
        args.min_query_coverage,
        args.max_subject_ratio
    )
    
    print(f"  Found {len(candidates)} frameshift candidates")
    
    # Write detailed text output
    with open(args.output, 'w') as out:
        out.write(f"FRAMESHIFT DETECTION RESULTS\n")
        out.write(f"{'='*80}\n\n")
        out.write(f"Parameters:\n")
        out.write(f"  Min identity: {args.min_identity}%\n")
        out.write(f"  Max genomic gap: {args.max_gap} bp\n")
        out.write(f"  Min query coverage: {args.min_query_coverage}%\n")
        out.write(f"  Max subject/query ratio: {args.max_subject_ratio}\n\n")
        out.write(f"Total candidates detected: {len(candidates)}\n")
        out.write(f"{'='*80}\n\n")
        
        for i, cand in enumerate(candidates, 1):
            out.write(f"{'*'*80}\n")
            out.write(f"CANDIDATE #{i}: POSSIBLE FRAMESHIFT DETECTED\n")
            out.write(f"{'*'*80}\n")
            out.write(f"Reference gene:   {cand['q_contig']}:{cand['q_start']}-{cand['q_end']} ({cand['q_len']} bp, strand {cand['q_strand']})\n")
            out.write(f"\n")
            out.write(f"Fragment 1:       {cand['s1_contig']}:{cand['s1_start']}-{cand['s1_end']} ({cand['s1_len']} bp)\n")
            out.write(f"  Query coverage: {cand['s1_qstart']}-{cand['s1_qend']} ({cand['s1_qend']-cand['s1_qstart']} bp)\n")
            out.write(f"  Identity:       {cand['s1_identity']:.3f}%\n")
            out.write(f"  Size ratio:     {cand['s1_len']/cand['q_len']:.2f}x query\n")
            out.write(f"\n")
            out.write(f"Fragment 2:       {cand['s2_contig']}:{cand['s2_start']}-{cand['s2_end']} ({cand['s2_len']} bp)\n")
            out.write(f"  Query coverage: {cand['s2_qstart']}-{cand['s2_qend']} ({cand['s2_qend']-cand['s2_qstart']} bp)\n")
            out.write(f"  Identity:       {cand['s2_identity']:.3f}%\n")
            out.write(f"  Size ratio:     {cand['s2_len']/cand['q_len']:.2f}x query\n")
            out.write(f"\n")
            out.write(f"Genomic gap:      {cand['genomic_gap']} bp")
            if cand['genomic_gap'] < 0:
                out.write(f" (overlapping)")
            out.write(f"\n")
            out.write(f"Query coverage:   {cand['query_coverage']} bp ({cand['query_coverage_pct']:.1f}%)\n")
            out.write(f"Query gap:        {cand['query_gap']} bp\n")
            out.write(f"Query overlap:    {cand['query_overlap']} bp\n")
            out.write(f"Combined length:  {cand['combined_len']} bp ({cand['combined_len']/cand['q_len']:.2f}x query)\n")
            out.write(f"Strand:           {cand['strand']}\n")
            out.write(f"\n")
    
    # Write TSV output
    with open(args.output_tsv, 'w') as out:
        out.write("query\tq_contig\tq_start\tq_end\tq_len\tq_strand\t"
                 "subject1\ts1_contig\ts1_start\ts1_end\ts1_len\ts1_qstart\ts1_qend\ts1_identity\ts1_ratio\t"
                 "subject2\ts2_contig\ts2_start\ts2_end\ts2_len\ts2_qstart\ts2_qend\ts2_identity\ts2_ratio\t"
                 "genomic_gap\tquery_coverage\tquery_coverage_pct\tquery_gap\tquery_overlap\t"
                 "combined_len\tcombined_ratio\tstrand\n")
        
        for cand in candidates:
            s1_ratio = cand['s1_len'] / cand['q_len']
            s2_ratio = cand['s2_len'] / cand['q_len']
            combined_ratio = cand['combined_len'] / cand['q_len']
            
            out.write(f"{cand['query']}\t{cand['q_contig']}\t{cand['q_start']}\t{cand['q_end']}\t{cand['q_len']}\t{cand['q_strand']}\t"
                     f"{cand['subject1']}\t{cand['s1_contig']}\t{cand['s1_start']}\t{cand['s1_end']}\t{cand['s1_len']}\t"
                     f"{cand['s1_qstart']}\t{cand['s1_qend']}\t{cand['s1_identity']:.3f}\t{s1_ratio:.3f}\t"
                     f"{cand['subject2']}\t{cand['s2_contig']}\t{cand['s2_start']}\t{cand['s2_end']}\t{cand['s2_len']}\t"
                     f"{cand['s2_qstart']}\t{cand['s2_qend']}\t{cand['s2_identity']:.3f}\t{s2_ratio:.3f}\t"
                     f"{cand['genomic_gap']}\t{cand['query_coverage']}\t{cand['query_coverage_pct']:.2f}\t"
                     f"{cand['query_gap']}\t{cand['query_overlap']}\t{cand['combined_len']}\t{combined_ratio:.3f}\t{cand['strand']}\n")
    
    print(f"\nResults saved:")
    print(f"  Detailed: {args.output}")
    print(f"  TSV: {args.output_tsv}")
    
    # Summary statistics
    if candidates:
        print(f"\n=== SUMMARY STATISTICS ===")
        
        # Genomic gap distribution
        genomic_gaps = [c['genomic_gap'] for c in candidates]
        print(f"\nGenomic gap sizes:")
        print(f"  Min: {min(genomic_gaps)} bp")
        print(f"  Max: {max(genomic_gaps)} bp")
        print(f"  Average: {sum(genomic_gaps)/len(genomic_gaps):.1f} bp")
        overlapping = sum(1 for g in genomic_gaps if g < 0)
        print(f"  Overlapping fragments: {overlapping}")
        
        # Query coverage distribution
        coverages = [c['query_coverage_pct'] for c in candidates]
        print(f"\nQuery coverage:")
        print(f"  Min: {min(coverages):.1f}%")
        print(f"  Max: {max(coverages):.1f}%")
        print(f"  Average: {sum(coverages)/len(coverages):.1f}%")
        
        # Query gap distribution
        query_gaps = [c['query_gap'] for c in candidates]
        print(f"\nQuery gaps (uncovered regions):")
        print(f"  Min: {min(query_gaps)} bp")
        print(f"  Max: {max(query_gaps)} bp")
        print(f"  Average: {sum(query_gaps)/len(query_gaps):.1f} bp")
        no_gap = sum(1 for g in query_gaps if g == 0)
        print(f"  Complete coverage (no gap): {no_gap}")
        
        # Identity distribution
        identities = [c['s1_identity'] for c in candidates] + [c['s2_identity'] for c in candidates]
        print(f"\nFragment identities:")
        print(f"  Min: {min(identities):.2f}%")
        print(f"  Max: {max(identities):.2f}%")
        print(f"  Average: {sum(identities)/len(identities):.2f}%")
        
        # Size ratios
        ratios = [(c['s1_len']/c['q_len']) for c in candidates] + [(c['s2_len']/c['q_len']) for c in candidates]
        print(f"\nFragment/Query size ratios:")
        print(f"  Min: {min(ratios):.2f}x")
        print(f"  Max: {max(ratios):.2f}x")
        print(f"  Average: {sum(ratios)/len(ratios):.2f}x")

        # Combined ratios
        combined_ratios = [c['combined_len']/c['q_len'] for c in candidates]
        print(f"\nCombined/Query size ratios:")
        print(f"  Min: {min(combined_ratios):.2f}x")
        print(f"  Max: {max(combined_ratios):.2f}x")
        print(f"  Average: {sum(combined_ratios)/len(combined_ratios):.2f}x")
        
        # Strand distribution
        plus_strand = sum(1 for c in candidates if c['strand'] == '+')
        minus_strand = sum(1 for c in candidates if c['strand'] == '-')
        print(f"\nBy strand:")
        print(f"  Plus (+): {plus_strand}")
        print(f"  Minus (-): {minus_strand}")
        
        # Contig distribution
        contigs = set(c['s1_contig'] for c in candidates)
        print(f"\nContigs involved: {len(contigs)}")
        
        # Show top 5 queries with most fragments
        from collections import Counter
        query_counts = Counter(c['query'] for c in candidates)
        if query_counts:
            print(f"\nTop queries with multiple frameshift candidates:")
            for query, count in query_counts.most_common(5):
                print(f"  {query}: {count} candidate(s)")
        
        # Fragment size distribution
        frag_sizes = [c['s1_len'] for c in candidates] + [c['s2_len'] for c in candidates]
        print(f"\nFragment sizes:")
        print(f"  Min: {min(frag_sizes)} bp")
        print(f"  Max: {max(frag_sizes)} bp")
        print(f"  Average: {sum(frag_sizes)/len(frag_sizes):.1f} bp")
        
        # Query size distribution
        query_sizes = [c['q_len'] for c in candidates]
        print(f"\nQuery sizes:")
        print(f"  Min: {min(query_sizes)} bp")
        print(f"  Max: {max(query_sizes)} bp")
        print(f"  Average: {sum(query_sizes)/len(query_sizes):.1f} bp")
        
    else:
        print("\nNo frameshift candidates detected with current parameters.")
        print("Consider adjusting:")
        print("  --min-identity (lower to be more permissive)")
        print("  --max-gap (higher to allow larger gaps)")
        print("  --min-query-coverage (lower to be more permissive)")
        print("  --max-subject-ratio (higher to allow larger fragments)")

if __name__ == "__main__":
    main()
