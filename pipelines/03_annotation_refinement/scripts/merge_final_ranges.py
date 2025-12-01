#!/usr/bin/env python3

"""
Merge ranges from multiple sources into a single sorted file.
Enriches missing strand, n_tools, and tools information from annotation results.
"""

import sys
import argparse
import re
import pandas as pd

def parse_hhr_filename(filename):
    """
    Parse HHR filename to extract start, end, strand, frame.
    Format: start-end_(strand)_frame.hhr
    Example: 125806-127659_(-)_4.hhr
    Returns: (start, end, strand, frame) or None if parsing fails
    """
    match = re.match(r'(\d+)-(\d+)_\(([+-])\)_(\d+)\.hhr', filename)
    if match:
        start = int(match.group(1))
        end = int(match.group(2))
        strand = match.group(3)
        frame = int(match.group(4))
        return (start, end, strand, frame)
    return None

def load_annotation_lookup(annotation_file):
    """
    Load final_overlap_results_combined.tsv and create lookup dictionary.
    Key: (start, end) tuple
    Value: dict with strand, n_tools, tools (NO frame)
    """
    df = pd.read_csv(annotation_file, sep='\t')
    lookup = {}
    
    for _, row in df.iterrows():
        # Range can be in columns DF1 or DF2
        
        # First range (DF1 columns)
        start1 = row['Start_DF1']
        end1 = row['End_DF1']
        strand1 = row['Strand_DF1']
        tools1 = row['Tools_DF1']
        
        # Skip if NA
        if pd.notna(start1) and pd.notna(end1):
            start1 = int(start1)
            end1 = int(end1)
            
            # Count number of tools (comma-separated)
            if pd.notna(tools1):
                n_tools1 = len([t.strip() for t in str(tools1).split(',')])
            else:
                n_tools1 = 0
            
            key1 = (start1, end1)
            if key1 not in lookup:
                lookup[key1] = {
                    'strand': strand1 if pd.notna(strand1) else 'unknown',
                    'n_tools': n_tools1,
                    'tools': tools1 if pd.notna(tools1) else 'unknown'
                }
        
        # Second range (DF2 columns) - if exists
        start2 = row['Start_DF2']
        end2 = row['End_DF2']
        
        if pd.notna(start2) and pd.notna(end2):
            start2 = int(start2)
            end2 = int(end2)
            strand2 = row['Strand_DF2']
            tools2 = row['Tools_DF2']
            
            # Count number of tools
            if pd.notna(tools2):
                n_tools2 = len([t.strip() for t in str(tools2).split(',')])
            else:
                n_tools2 = 0
            
            key2 = (start2, end2)
            if key2 not in lookup:
                lookup[key2] = {
                    'strand': strand2 if pd.notna(strand2) else 'unknown',
                    'n_tools': n_tools2,
                    'tools': tools2 if pd.notna(tools2) else 'unknown'
                }
    
    return lookup

def enrich_range_info(ranges, annotation_lookup):
    """
    Enrich ranges with strand, n_tools, and tools from annotation lookup.
    Frame is kept as-is (unknown if not available).
    """
    enriched = []
    not_found = 0
    
    for r in ranges:
        key = (r['start'], r['end'])
        
        # Try to get info from annotation lookup
        if key in annotation_lookup:
            info = annotation_lookup[key]
            r['strand'] = info['strand']
            r['n_tools'] = info['n_tools']
            r['tools'] = info['tools']
            # Frame stays as-is (from original source or 'unknown')
        else:
            # Keep existing values or set defaults
            if r['strand'] == 'unknown':
                not_found += 1
            if 'n_tools' not in r or r['n_tools'] == 'unknown':
                r['n_tools'] = 'unknown'
            if 'tools' not in r or r['tools'] == 'unknown':
                r['tools'] = 'unknown'
        
        enriched.append(r)
    
    if not_found > 0:
        print(f"  Warning: {not_found} ranges not found in annotation file", file=sys.stderr)
    
    return enriched

def load_cat2_ranges(cat2_file):
    """
    Load ranges from cat2_full_consensus_5.tsv
    Columns: Start, End, n_tools_observed, tools_observed
    """
    df = pd.read_csv(cat2_file, sep='\t')
    ranges = []
    
    for _, row in df.iterrows():
        ranges.append({
            'start': int(row['Start']),
            'end': int(row['End']),
            'strand': 'unknown',
            'frame': 'unknown',
            'source': 'cat2_consensus_5',
            'n_tools': row['n_tools_observed'],
            'tools': row['tools_observed']
        })
    
    return ranges

def load_best_query_ranges(best_query_file):
    """
    Load ranges from best_query_per_group.tsv
    Parse the 'best_query' column (format: start-end_(strand)_frame.hhr)
    """
    df = pd.read_csv(best_query_file, sep='\t')
    ranges = []
    
    for _, row in df.iterrows():
        best_query = row['best_query']
        parsed = parse_hhr_filename(best_query)
        
        if parsed:
            start, end, strand, frame = parsed
            ranges.append({
                'start': start,
                'end': end,
                'strand': strand,
                'frame': frame,
                'source': 'best_query_group',
                'n_tools': 'unknown',  # Will be enriched later
                'tools': 'unknown'
            })
        else:
            print(f"Warning: Could not parse best_query: {best_query}", file=sys.stderr)
    
    return ranges

def load_filtered_ranges(filtered_file):
    """
    Load ranges from all_ranges_to_keep.tsv
    Parse the 'file' column (format: start-end_(strand)_frame.hhr)
    """
    df = pd.read_csv(filtered_file, sep='\t')
    ranges = []
    
    for _, row in df.iterrows():
        filename = row['file']
        parsed = parse_hhr_filename(filename)
        
        if parsed:
            start, end, strand, frame = parsed
            ranges.append({
                'start': start,
                'end': end,
                'strand': strand,
                'frame': frame,
                'source': 'filtered_singles',
                'n_tools': row['n_annotation_tools'],
                'tools': row['annotation_tools']
            })
        else:
            print(f"Warning: Could not parse file: {filename}", file=sys.stderr)
    
    return ranges

def remove_duplicates(ranges):
    """
    Remove duplicate ranges based on (start, end).
    Keep the first occurrence with most complete information.
    """
    seen = {}
    unique_ranges = []
    
    for r in ranges:
        key = (r['start'], r['end'])
        
        if key not in seen:
            seen[key] = r
            unique_ranges.append(r)
        else:
            # If new range has better info (not unknown), update
            existing = seen[key]
            
            # Update strand if existing is unknown
            if existing['strand'] == 'unknown' and r['strand'] != 'unknown':
                existing['strand'] = r['strand']
            
            # Update frame if existing is unknown
            if existing['frame'] == 'unknown' and r['frame'] != 'unknown':
                existing['frame'] = r['frame']
            
            # Update n_tools if existing is unknown
            if existing['n_tools'] == 'unknown' and r['n_tools'] != 'unknown':
                existing['n_tools'] = r['n_tools']
            
            # Update tools if existing is unknown
            if existing['tools'] == 'unknown' and r['tools'] != 'unknown':
                existing['tools'] = r['tools']
            
            # Track multiple sources
            if 'sources' not in existing:
                existing['sources'] = [existing['source']]
            if r['source'] not in existing['sources']:
                existing['sources'].append(r['source'])
    
    # Update source field for duplicates
    for r in unique_ranges:
        if 'sources' in r:
            r['source'] = ';'.join(r['sources'])
    
    return unique_ranges

def parse_args():
    parser = argparse.ArgumentParser(
        description='Merge ranges from multiple annotation sources'
    )
    parser.add_argument('--cat2', required=True, 
                       help='cat2_full_consensus_5.tsv file')
    parser.add_argument('--best-query', required=True,
                       help='best_query_per_group.tsv file')
    parser.add_argument('--filtered', required=True,
                       help='all_ranges_to_keep.tsv file')
    parser.add_argument('--annotation', required=True,
                       help='final_overlap_results_combined.tsv file')
    parser.add_argument('--chrom-name', default='Lokiarchaeum_B36',
                       help='Chromosome/contig name for BED file (default: Lokiarchaeum_B36)')
    parser.add_argument('--output', default='merged_final_ranges.tsv',
                       help='Output file (default: merged_final_ranges.tsv)')
    parser.add_argument('--bed-output', default='merged_final_ranges.bed',
                       help='Output BED file (default: merged_final_ranges.bed)')
    return parser.parse_args()

def main():
    args = parse_args()
    
    print("Loading annotation lookup table...")
    annotation_lookup = load_annotation_lookup(args.annotation)
    print(f"  Loaded {len(annotation_lookup)} unique ranges from annotation file")
    
    print("\nLoading ranges from cat2_full_consensus_5.tsv...")
    cat2_ranges = load_cat2_ranges(args.cat2)
    print(f"  Loaded {len(cat2_ranges)} ranges")
    
    print("Loading ranges from best_query_per_group.tsv...")
    best_query_ranges = load_best_query_ranges(args.best_query)
    print(f"  Loaded {len(best_query_ranges)} ranges")
    
    print("Loading ranges from all_ranges_to_keep.tsv...")
    filtered_ranges = load_filtered_ranges(args.filtered)
    print(f"  Loaded {len(filtered_ranges)} ranges")
    
    # Combine all ranges
    all_ranges = cat2_ranges + best_query_ranges + filtered_ranges
    print(f"\nTotal ranges before deduplication: {len(all_ranges)}")
    
    # Remove duplicates
    unique_ranges = remove_duplicates(all_ranges)
    print(f"Unique ranges after deduplication: {len(unique_ranges)}")
    
    # Enrich with annotation information (strand, n_tools, tools - NOT frame)
    print("\nEnriching ranges with annotation information...")
    enriched_ranges = enrich_range_info(unique_ranges, annotation_lookup)
    
    # Count how many were enriched
    enriched_strand = sum(1 for r in enriched_ranges if r['strand'] != 'unknown')
    enriched_tools = sum(1 for r in enriched_ranges if r['n_tools'] != 'unknown')
    print(f"  Ranges with strand information: {enriched_strand}/{len(enriched_ranges)}")
    print(f"  Ranges with tools information: {enriched_tools}/{len(enriched_ranges)}")
    
    # Sort by start position
    enriched_ranges.sort(key=lambda x: x['start'])
    
    # Count how many were enriched
    enriched_strand = sum(1 for r in enriched_ranges if r['strand'] != 'unknown')
    enriched_tools = sum(1 for r in enriched_ranges if r['n_tools'] != 'unknown')
    print(f"  Ranges with strand information: {enriched_strand}/{len(enriched_ranges)}")
    print(f"  Ranges with tools information: {enriched_tools}/{len(enriched_ranges)}")
    
    # Sort by start position
    enriched_ranges.sort(key=lambda x: x['start'])
    
    # Create DataFrame with only desired columns
    df = pd.DataFrame(enriched_ranges)
    
    # Select and reorder columns: start, end, strand, frame, source, n_tools, tools
    output_columns = ['start', 'end', 'strand', 'frame', 'source', 'n_tools', 'tools']
    df_output = df[output_columns]
    
    # Save to TSV
    df_output.to_csv(args.output, sep='\t', index=False)
    print(f"\nMerged ranges saved to: {args.output}")
    
    # Create BED file (0-based, half-open intervals)
    print(f"\nCreating BED file with chromosome name: {args.chrom_name}")
    with open(args.bed_output, 'w') as bed:
        for _, row in df_output.iterrows():
            # BED format: chrom, chromStart (0-based), chromEnd, name, score, strand
            chrom = args.chrom_name
            start_0based = row['start'] - 1
            end_pos = row['end']
            
            # Create name: start-end_strand_frame
            strand_str = row['strand'] if row['strand'] != 'unknown' else '.'
            frame_str = str(row['frame']) if row['frame'] != 'unknown' else 'unk'
            name = f"{row['start']}-{row['end']}_{strand_str}_{frame_str}"
            
            score = "0"
            bed_strand = row['strand'] if row['strand'] in ['+', '-'] else '.'
            
            bed.write(f"{chrom}\t{start_0based}\t{end_pos}\t{name}\t{score}\t{bed_strand}\n")
    
    print(f"BED file saved to: {args.bed_output}")
    
    # Summary statistics
    print("\n=== SUMMARY ===")
    print(f"Total unique ranges: {len(df_output)}")
    
    print("\nRanges by source:")
    print(df_output['source'].value_counts().to_string())
    
    print("\nRanges by strand:")
    strand_counts = df_output['strand'].value_counts()
    for strand, count in strand_counts.items():
        pct = (count / len(df_output)) * 100
        print(f"  {strand}: {count} ({pct:.1f}%)")
    
    print("\nRanges by n_tools:")
    tools_counts = df_output['n_tools'].value_counts().sort_index()
    for n_tools, count in tools_counts.items():
        pct = (count / len(df_output)) * 100
        print(f"  {n_tools} tools: {count} ({pct:.1f}%)")
    
    # Genome coverage statistics
    total_length = (df_output['end'] - df_output['start']).sum()
    genome_span = df_output['end'].max() - df_output['start'].min() + 1
    
    print(f"\nGenome span (first to last range): {genome_span:,} bp")
    print(f"Total bases covered by ranges: {total_length:,} bp")
    print(f"Coding density: {(total_length/genome_span)*100:.1f}%")
    
    # Frame statistics (excluding unknown)
    frame_known = df_output[df_output['frame'] != 'unknown']
    if len(frame_known) > 0:
        print(f"\nRanges with known frame: {len(frame_known)}/{len(df_output)}")
        print("Frame distribution (known only):")
        print(frame_known['frame'].value_counts().sort_index().to_string())
    else:
        print(f"\nRanges with known frame: 0/{len(df_output)}")

if __name__ == "__main__":
    main()
