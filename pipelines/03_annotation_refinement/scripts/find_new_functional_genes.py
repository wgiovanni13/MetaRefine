#!/usr/bin/env python3
"""
Find new functional genes not predicted by a specific tool, not overlapping with that tool,
and having functional annotation in UniRef90
"""

import sys
import gzip
from pathlib import Path
from collections import defaultdict

def is_characterized(description):
    """Check if protein is characterized"""
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

def extract_tool_ranges_from_merged(merged_file, tool_name):
    """Extract all ranges predicted by a specific tool from merged table"""
    tool_ranges = []
    
    with open(merged_file, 'r') as f:
        header = f.readline()
        
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            
            start_df1 = fields[0]
            end_df1 = fields[1]
            start_df2 = fields[2]
            end_df2 = fields[3]
            tools_df1 = fields[10]
            tools_df2 = fields[11]
            
            # Check DF1
            if start_df1 != "NA" and end_df1 != "NA":
                if tool_name in tools_df1.lower():
                    tool_ranges.append((int(start_df1), int(end_df1)))
            
            # Check DF2
            if start_df2 != "NA" and end_df2 != "NA":
                if tool_name in tools_df2.lower():
                    tool_ranges.append((int(start_df2), int(end_df2)))
    
    # Remove duplicates
    tool_ranges = list(set(tool_ranges))
    tool_ranges.sort()
    
    print(f"[INFO] Extracted {len(tool_ranges)} {tool_name.upper()} gene ranges", file=sys.stderr)
    return tool_ranges

def has_overlap_with_tool(start, end, tool_ranges):
    """Check if a range overlaps with any tool range"""
    for t_start, t_end in tool_ranges:
        # Check for any overlap
        if not (end < t_start or start > t_end):
            return True
    return False

def parse_merged_table(merged_file, tool_ranges, tool_name):
    """Parse merged annotation table and extract ranges NOT predicted by tool and NOT overlapping"""
    non_tool_ranges = []
    
    with open(merged_file, 'r') as f:
        header = f.readline()
        
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            
            start_df1 = fields[0]
            end_df1 = fields[1]
            start_df2 = fields[2]
            end_df2 = fields[3]
            strand_df1 = fields[6]
            strand_df2 = fields[7]
            tools_df1 = fields[10]
            tools_df2 = fields[11]
            
            # Check DF1 (columns 1-2)
            if start_df1 != "NA" and end_df1 != "NA":
                if tool_name not in tools_df1.lower():
                    start = int(start_df1)
                    end = int(end_df1)
                    
                    # Check if it overlaps with any tool range
                    if not has_overlap_with_tool(start, end, tool_ranges):
                        non_tool_ranges.append({
                            'start': start,
                            'end': end,
                            'strand': strand_df1,
                            'tools': tools_df1,
                            'source': 'DF1'
                        })
            
            # Check DF2 (columns 3-4)
            if start_df2 != "NA" and end_df2 != "NA":
                if tool_name not in tools_df2.lower():
                    start = int(start_df2)
                    end = int(end_df2)
                    
                    # Check if it overlaps with any tool range
                    if not has_overlap_with_tool(start, end, tool_ranges):
                        non_tool_ranges.append({
                            'start': start,
                            'end': end,
                            'strand': strand_df2,
                            'tools': tools_df2,
                            'source': 'DF2'
                        })
    
    print(f"[INFO] Found {len(non_tool_ranges)} ranges NOT predicted by {tool_name.upper()} and NOT overlapping", file=sys.stderr)
    
    # Remove duplicates and contained ranges
    non_tool_ranges = remove_contained_ranges(non_tool_ranges)
    print(f"[INFO] After removing duplicates/contained: {len(non_tool_ranges)} ranges", file=sys.stderr)
    
    return non_tool_ranges

def remove_contained_ranges(ranges):
    """Remove duplicate ranges and ranges contained within others"""
    if not ranges:
        return []
    
    # Sort by start, then by length (descending)
    ranges.sort(key=lambda x: (x['start'], -(x['end'] - x['start'])))
    
    # Remove exact duplicates
    unique_ranges = []
    seen = set()
    
    for r in ranges:
        key = (r['start'], r['end'])
        if key not in seen:
            seen.add(key)
            unique_ranges.append(r)
    
    # Remove contained ranges
    filtered_ranges = []
    
    for i, current in enumerate(unique_ranges):
        is_contained = False
        
        for other in unique_ranges:
            if current is other:
                continue
            
            # Check if current is contained in other
            if (current['start'] >= other['start'] and 
                current['end'] <= other['end'] and
                (current['start'] != other['start'] or current['end'] != other['end'])):
                is_contained = True
                break
        
        if not is_contained:
            filtered_ranges.append(current)
    
    return filtered_ranges

def parse_diamond_chunk(diamond_file):
    """Parse DIAMOND results and extract characterized hits"""
    hits = defaultdict(list)
    
    opener = gzip.open if diamond_file.endswith('.gz') else open
    
    with opener(diamond_file, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 15:
                continue
            
            qseqid = fields[0]
            sseqid = fields[1]
            pident = float(fields[2])
            length = int(fields[3])
            evalue = float(fields[10])
            bitscore = float(fields[11])
            qlen = int(fields[12])
            slen = int(fields[13])
            stitle = fields[14] if len(fields) > 14 else "No description"
            
            # Extract coordinates from query ID
            try:
                coords = qseqid.split('(')[0]
                start, end = coords.split('-')
                start = int(start)
                end = int(end)
                
                # Extract strand
                strand = qseqid.split('(')[1].split(')')[0]
                
                query_cov = (length / qlen) * 100
                
                # Only keep characterized proteins with good hits
                if is_characterized(stitle) and evalue < 1e-5 and query_cov >= 50:
                    hits[(start, end)].append({
                        'ref_id': sseqid,
                        'description': stitle,
                        'identity': pident,
                        'coverage': query_cov,
                        'evalue': evalue,
                        'bitscore': bitscore,
                        'strand': strand
                    })
            except (ValueError, IndexError):
                continue
    
    return hits

def find_new_functional_genes(merged_file, diamond_dir, output_file, tool_name):
    """Main function to find new functional genes"""
    
    # 1. Extract tool ranges from merged table
    print(f"[INFO] Extracting {tool_name.upper()} ranges from merged table...", file=sys.stderr)
    tool_ranges = extract_tool_ranges_from_merged(merged_file, tool_name)
    
    # 2. Parse merged table for non-tool ranges
    print("[INFO] Parsing merged annotation table...", file=sys.stderr)
    non_tool_ranges = parse_merged_table(merged_file, tool_ranges, tool_name)
    
    if len(non_tool_ranges) == 0:
        print(f"[WARNING] No non-{tool_name.upper()}, non-overlapping ranges found!", file=sys.stderr)
        # Create empty output file
        with open(output_file, 'w') as out:
            out.write('\t'.join([
                'gene_range', 'strand', 'length_bp', 'predicted_by_tools',
                'ref_protein_id', 'protein_description', 'identity_%',
                'coverage_%', 'evalue', 'bitscore'
            ]) + '\n')
        return
    
    # 3. Parse all DIAMOND chunks
    print("[INFO] Parsing DIAMOND results...", file=sys.stderr)
    all_hits = defaultdict(list)
    
    diamond_files = list(Path(diamond_dir).glob("*.diamond.gz"))
    print(f"[INFO] Found {len(diamond_files)} DIAMOND chunk files", file=sys.stderr)
    
    for diamond_file in diamond_files:
        chunk_hits = parse_diamond_chunk(str(diamond_file))
        for coords, hits in chunk_hits.items():
            all_hits[coords].extend(hits)
    
    print(f"[INFO] Found hits for {len(all_hits)} coordinate ranges", file=sys.stderr)
    
    # 4. Match non-tool ranges with DIAMOND hits
    print("[INFO] Matching ranges with functional annotations...", file=sys.stderr)
    new_functional_genes = []
    
    for range_info in non_tool_ranges:
        start = range_info['start']
        end = range_info['end']
        
        # Find matching hits (exact or overlapping coordinates)
        matching_hits = []
        for (hit_start, hit_end), hits in all_hits.items():
            # Check if ranges overlap significantly
            overlap_start = max(start, hit_start)
            overlap_end = min(end, hit_end)
            overlap_len = max(0, overlap_end - overlap_start + 1)
            
            range_len = end - start + 1
            overlap_pct = (overlap_len / range_len) * 100
            
            if overlap_pct >= 80:  # At least 80% overlap
                matching_hits.extend(hits)
        
        if matching_hits:
            # Sort by bitscore and take best hit
            matching_hits.sort(key=lambda x: x['bitscore'], reverse=True)
            best_hit = matching_hits[0]
            
            new_functional_genes.append({
                'range': f"{start}-{end}",
                'strand': range_info['strand'],
                'tools': range_info['tools'],
                'length': end - start + 1,
                'ref_id': best_hit['ref_id'],
                'description': best_hit['description'],
                'identity': best_hit['identity'],
                'coverage': best_hit['coverage'],
                'evalue': best_hit['evalue'],
                'bitscore': best_hit['bitscore']
            })
    
    # 5. Write output
    print(f"[INFO] Found {len(new_functional_genes)} new functional genes", file=sys.stderr)
    
    with open(output_file, 'w') as out:
        # Header
        out.write('\t'.join([
            'gene_range',
            'strand',
            'length_bp',
            'predicted_by_tools',
            'ref_protein_id',
            'protein_description',
            'identity_%',
            'coverage_%',
            'evalue',
            'bitscore'
        ]) + '\n')
        
        # Sort by bitscore
        new_functional_genes.sort(key=lambda x: x['bitscore'], reverse=True)
        
        for gene in new_functional_genes:
            out.write('\t'.join([
                gene['range'],
                gene['strand'],
                str(gene['length']),
                gene['tools'],
                gene['ref_id'],
                gene['description'],
                f"{gene['identity']:.1f}",
                f"{gene['coverage']:.1f}",
                f"{gene['evalue']:.2e}",
                f"{gene['bitscore']:.1f}"
            ]) + '\n')
    
    print(f"[INFO] Results written to {output_file}", file=sys.stderr)
    
    # Print summary
    print("\n" + "=" * 80, file=sys.stderr)
    print(f"SUMMARY: New Functional Genes Not in {tool_name.upper()}", file=sys.stderr)
    print("=" * 80, file=sys.stderr)
    print(f"{tool_name.upper()} ranges: {len(tool_ranges)}", file=sys.stderr)
    print(f"Non-{tool_name.upper()}, non-overlapping ranges: {len(non_tool_ranges)}", file=sys.stderr)
    print(f"New functional genes found: {len(new_functional_genes)}", file=sys.stderr)
    
    if new_functional_genes:
        print(f"\nTop 10 new functional genes (by bitscore):", file=sys.stderr)
        for i, gene in enumerate(new_functional_genes[:10], 1):
            desc_short = gene['description'][:55] + "..." if len(gene['description']) > 55 else gene['description']
            print(f"\n{i}. Range: {gene['range']} ({gene['length']}bp) | Strand: {gene['strand']}", file=sys.stderr)
            print(f"   Description: {desc_short}", file=sys.stderr)
            print(f"   Identity: {gene['identity']:.1f}% | Coverage: {gene['coverage']:.1f}%", file=sys.stderr)
            print(f"   Predicted by: {gene['tools']}", file=sys.stderr)
    
    print("=" * 80, file=sys.stderr)

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("=" * 80, file=sys.stderr)
        print("Find New Functional Genes", file=sys.stderr)
        print("=" * 80, file=sys.stderr)
        print("\nUsage: find_new_functional_genes.py <merged_table.tsv> <diamond_results_dir> <output.tsv> <tool_name>", file=sys.stderr)
        print("\nExamples:", file=sys.stderr)
        print("  # Find genes not in Prokka:", file=sys.stderr)
        print("  python3 find_new_functional_genes.py merged.tsv diamond_dir/ new_genes_prokka.tsv prokka", file=sys.stderr)
        print("\n  # Find genes not in Prodigal:", file=sys.stderr)
        print("  python3 find_new_functional_genes.py merged.tsv diamond_dir/ new_genes_prodigal.tsv prodigal", file=sys.stderr)
        print("\n  # Find genes not in GMS2:", file=sys.stderr)
        print("  python3 find_new_functional_genes.py merged.tsv diamond_dir/ new_genes_gms2.tsv gms2", file=sys.stderr)
        print("\nAvailable tools: prokka, prodigal, gms2, glimmer, rast", file=sys.stderr)
        print("=" * 80, file=sys.stderr)
        sys.exit(1)
    
    merged_file = sys.argv[1]
    diamond_dir = sys.argv[2]
    output_file = sys.argv[3]
    tool_name = sys.argv[4].lower()
    
    find_new_functional_genes(merged_file, diamond_dir, output_file, tool_name)

