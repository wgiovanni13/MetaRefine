#!/usr/bin/env python3

"""
Analyze singles.tsv to find overlaps and annotation tool counts.
"""

import sys
import argparse
import pandas as pd

def calculate_overlap(start1, end1, start2, end2):

    """
    Calculate overlap between two ranges.
    Returns: (overlap_length, overlap_percentage_of_first, overlap_percentage_of_second)
    """

    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    
    if overlap_start >= overlap_end:
        return (0, 0.0, 0.0)
    
    overlap_length = overlap_end - overlap_start
    len1 = end1 - start1
    len2 = end2 - start2
    
    pct1 = (overlap_length / len1 * 100) if len1 > 0 else 0
    pct2 = (overlap_length / len2 * 100) if len2 > 0 else 0
    
    return (overlap_length, pct1, pct2)

def find_annotation_tools(start, end, strand, annotation_df):

    """
    Find how many tools predicted a given range.
    Returns: (n_tools, tools_list)
    """

    if annotation_df is None:
        return (0, "annotation_file_not_found")
    
    # Search in DF1 columns (Start_DF1, End_DF1)

    mask_df1 = (
        (annotation_df['Start_DF1'] == start) & 
        (annotation_df['End_DF1'] == end) &
        (annotation_df['Strand_DF1'] == strand)
    )
    
    # Search in DF2 columns (Start_DF2, End_DF2)

    mask_df2 = (
        (annotation_df['Start_DF2'] == start) & 
        (annotation_df['End_DF2'] == end) &
        (annotation_df['Strand_DF2'] == strand)
    )
    
    # Check DF1 first

    if mask_df1.any():
        row = annotation_df[mask_df1].iloc[0]
        tools = row['Tools_DF1']
        if pd.notna(tools):
            tools_list = [t.strip() for t in str(tools).split(',')]
            return (len(tools_list), tools)
    
    # Check DF2

    if mask_df2.any():
        row = annotation_df[mask_df2].iloc[0]
        tools = row['Tools_DF2']
        if pd.notna(tools):
            tools_list = [t.strip() for t in str(tools).split(',')]
            return (len(tools_list), tools)
    
    return (0, "not_found_in_annotations")

def parse_args():

    """Parse command line arguments."""

    parser = argparse.ArgumentParser(
        description='Analyze singles.tsv for overlaps and annotation tool counts'
    )
    parser.add_argument(
        '--singles',
        required=True,
        help='Path to singles.tsv file'
    )
    parser.add_argument(
        '--annotation',
        required=True,
        help='Path to annotation TSV file (final_overlap_results_combined.tsv)'
    )
    parser.add_argument(
        '--outdir',
        default='.',
        help='Output directory (default: current directory)'
    )
    return parser.parse_args()

def main():

    args = parse_args()
    
    print("Loading singles.tsv...")
    try:
        singles_df = pd.read_csv(args.singles, sep='\t')
        print(f"Loaded {len(singles_df)} single ranges")
    except Exception as e:
        print(f"ERROR: Could not load singles file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print("Loading annotation data...")
    try:
        annotation_df = pd.read_csv(args.annotation, sep='\t')
        print(f"Loaded {len(annotation_df)} annotation records")
    except Exception as e:
        print(f"WARNING: Could not load annotation file: {e}", file=sys.stderr)
        annotation_df = None
    
    # Prepare results

    results = []
    
    print("Analyzing overlaps and annotation tools...")
    total = len(singles_df)
    
    for idx, row in singles_df.iterrows():
        if (idx + 1) % 100 == 0:
            print(f"Processing {idx + 1}/{total}...")
        
        file1 = row['file']
        start1 = row['start']
        end1 = row['end']
        strand1 = row['strand']
        frame1 = row['frame']
        
        # Get annotation tools for this range

        n_tools, tools_str = find_annotation_tools(start1, end1, strand1, annotation_df)
        
        # Find overlaps with other ranges

        overlapping_ranges = []
        overlap_details = []
        
        for idx2, row2 in singles_df.iterrows():
            if idx == idx2:
                continue
            
            start2 = row2['start']
            end2 = row2['end']
            strand2 = row2['strand']
            
            overlap_len, pct1, pct2 = calculate_overlap(start1, end1, start2, end2)
            
            if overlap_len > 0:
                overlapping_ranges.append(row2['file'])
                overlap_details.append(
                    f"{row2['file']}({strand2}:{start2}-{end2}|{pct1:.1f}%)"
                )
        
        # Add result

        results.append({
            'file': file1,
            'strand': strand1,
            'start': start1,
            'end': end1,
            'frame': frame1,
            'length': end1 - start1,
            'n_annotation_tools': n_tools,
            'annotation_tools': tools_str,
            'n_overlapping_ranges': len(overlapping_ranges),
            'overlapping_files': ';'.join(overlapping_ranges) if overlapping_ranges else 'none',
            'overlap_details': ';'.join(overlap_details) if overlap_details else 'none'
        })
    
    # Create output DataFrame

    results_df = pd.DataFrame(results)
    
    # Save detailed results

    output_file = f"{args.outdir}/singles_detailed_analysis.tsv"
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nDetailed analysis saved to: {output_file}")
    
    # Generate summary statistics

    summary_file = f"{args.outdir}/singles_analysis_summary.tsv"
    with open(summary_file, 'w') as f:
        f.write("=== SINGLES ANALYSIS SUMMARY ===\n\n")
        f.write(f"Total single ranges analyzed: {len(results_df)}\n\n")
        
        # Annotation tools distribution

        f.write("Distribution by number of annotation tools:\n")
        tool_counts = results_df['n_annotation_tools'].value_counts().sort_index()
        for n_tools, count in tool_counts.items():
            f.write(f"  {n_tools} tools: {count} ranges\n")
        
        f.write("\n")
        
        # Overlap distribution

        f.write("Distribution by number of overlapping ranges:\n")
        overlap_counts = results_df['n_overlapping_ranges'].value_counts().sort_index()
        for n_overlaps, count in overlap_counts.items():
            f.write(f"  {n_overlaps} overlaps: {count} ranges\n")
        
        f.write("\n")
        
        # Ranges with no annotation tools

        no_tools = results_df[results_df['n_annotation_tools'] == 0]
        f.write(f"Ranges with NO annotation tools: {len(no_tools)}\n")
        
        # Ranges with overlaps but no annotation

        overlaps_no_tools = results_df[
            (results_df['n_overlapping_ranges'] > 0) & 
            (results_df['n_annotation_tools'] == 0)
        ]
        f.write(f"Ranges with overlaps but NO annotation tools: {len(overlaps_no_tools)}\n")
        
        # Ranges with annotation but no overlaps

        tools_no_overlaps = results_df[
            (results_df['n_annotation_tools'] > 0) & 
            (results_df['n_overlapping_ranges'] == 0)
        ]
        f.write(f"Ranges with annotation tools but NO overlaps: {len(tools_no_overlaps)}\n")
        
        # Ranges with both annotation and overlaps

        both = results_df[
            (results_df['n_annotation_tools'] > 0) & 
            (results_df['n_overlapping_ranges'] > 0)
        ]
        f.write(f"Ranges with BOTH annotation tools AND overlaps: {len(both)}\n")
    
    print(f"Summary statistics saved to: {summary_file}")
    
    # Print summary to console

    print("\n=== SUMMARY ===")
    print(f"Total single ranges: {len(results_df)}")
    print(f"Ranges with annotation tools: {len(results_df[results_df['n_annotation_tools'] > 0])}")
    print(f"Ranges with overlaps: {len(results_df[results_df['n_overlapping_ranges'] > 0])}")
    print(f"Ranges with NO annotation: {len(results_df[results_df['n_annotation_tools'] == 0])}")

if __name__ == "__main__":
    main()

