#!/usr/bin/env python3

"""
Filter singles based on annotation tools, overlaps, and coverage by groups.
"""

import sys
import argparse
import pandas as pd
import re

def parse_group_ranges(groups_file):
    """
    Parse groups.tsv and extract all individual ranges.
    Returns: list of tuples (start, end, strand, frame)
    """
    df = pd.read_csv(groups_file, sep='\t')
    all_ranges = []
    
    for _, row in df.iterrows():
        files = row['files'].split(';')
        strand = row['strand']
        
        for file in files:
            # Parse: start-end_(strand)_frame.hhr
            match = re.match(r'(\d+)-(\d+)_\(([+-])\)_(\d+)\.hhr', file)
            if match:
                start = int(match.group(1))
                end = int(match.group(2))
                strand = match.group(3)
                frame = int(match.group(4))
                all_ranges.append((start, end, strand, frame, file))
    
    return all_ranges

def calculate_coverage_percentage(single_start, single_end, group_start, group_end):
    """
    Calculate what percentage of the single range is covered by the group range.
    """
    overlap_start = max(single_start, group_start)
    overlap_end = min(single_end, group_end)
    
    if overlap_start >= overlap_end:
        return 0.0
    
    overlap_length = overlap_end - overlap_start
    single_length = single_end - single_start
    
    return (overlap_length / single_length * 100) if single_length > 0 else 0.0

def is_covered_by_groups(single_start, single_end, group_ranges, threshold=90.0):
    """
    Check if a single range is covered >= threshold% by any group range.
    Returns: (is_covered, covering_file, coverage_percentage)
    """
    for g_start, g_end, g_strand, g_frame, g_file in group_ranges:
        coverage = calculate_coverage_percentage(single_start, single_end, g_start, g_end)
        if coverage >= threshold:
            return (True, g_file, coverage)
    return (False, None, 0.0)

def get_five_tool_ranges(detailed_file):
    """
    Get all ranges predicted by 5 tools.
    Returns: set of filenames
    """
    df = pd.read_csv(detailed_file, sep='\t')
    five_tools = df[df['n_annotation_tools'] == 5]['file'].tolist()
    return set(five_tools)

def overlaps_only_with_five_tools(row, five_tool_ranges):
    """
    Check if a range only overlaps with 5-tool ranges.
    """
    if row['n_overlapping_ranges'] == 0:
        return False
    
    if row['overlapping_files'] == 'none':
        return False
    
    overlapping = row['overlapping_files'].split(';')
    
    # Check if ALL overlapping files are in five_tool_ranges
    return all(f in five_tool_ranges for f in overlapping)

def has_different_strand_overlaps(row):
    """
    Check if range has overlaps with different strand.
    """
    if row['n_overlapping_ranges'] == 0:
        return False
    
    if row['overlap_details'] == 'none':
        return False
    
    my_strand = row['strand']
    details = row['overlap_details'].split(';')
    
    for detail in details:
        # Parse: filename(strand:start-end|pct%)
        match = re.search(r'\(([+-]):', detail)
        if match:
            other_strand = match.group(1)
            if other_strand != my_strand:
                return True
    
    return False

def has_same_strand_low_overlap(row, threshold=30.0):
    """
    Check if range has overlaps with SAME strand but < threshold% overlap.
    Returns True if ALL same-strand overlaps are < threshold%
    """
    if row['n_overlapping_ranges'] == 0:
        return False
    
    if row['overlap_details'] == 'none':
        return False
    
    my_strand = row['strand']
    details = row['overlap_details'].split(';')
    
    has_same_strand = False
    all_low_overlap = True
    
    for detail in details:
        # Parse: filename(strand:start-end|pct%)
        match = re.search(r'\(([+-]):(\d+)-(\d+)\|([0-9.]+)%\)', detail)
        if match:
            other_strand = match.group(1)
            overlap_pct = float(match.group(4))
            
            if other_strand == my_strand:
                has_same_strand = True
                if overlap_pct >= threshold:
                    all_low_overlap = False
                    break
    
    return has_same_strand and all_low_overlap

def has_no_overlaps(row):
    """
    Check if range has NO overlaps at all.
    """
    return row['n_overlapping_ranges'] == 0

def check_asymmetric_overlap(row, detailed_df, high_threshold=80.0, low_threshold=45.0):
    """
    Check if this range should be removed due to asymmetric overlap:
    - This range is covered >high_threshold% by another range
    - This range covers <=low_threshold% of that other range
    - Both have same strand
    
    Returns: (should_remove, other_file, my_coverage, other_coverage)
    """
    if row['n_overlapping_ranges'] == 0:
        return (False, None, 0.0, 0.0)
    
    if row['overlap_details'] == 'none':
        return (False, None, 0.0, 0.0)
    
    my_strand = row['strand']
    my_file = row['file']
    details = row['overlap_details'].split(';')
    
    for detail in details:
        # Parse: filename(strand:start-end|pct%)
        match = re.match(r'(.+?)\(([+-]):(\d+)-(\d+)\|([0-9.]+)%\)', detail.strip())
        if match:
            other_file = match.group(1)
            other_strand = match.group(2)
            my_coverage_pct = float(match.group(5))
            
            # Only check same strand overlaps
            if other_strand != my_strand:
                continue
            
            # Check if I'm covered >high_threshold% by the other range
            if my_coverage_pct > high_threshold:
                # Now check how much I cover of the other range
                other_row = detailed_df[detailed_df['file'] == other_file]
                
                if not other_row.empty:
                    other_row = other_row.iloc[0]
                    other_details = other_row['overlap_details']
                    
                    if other_details != 'none':
                        # Find my coverage in the other range's details
                        other_match = re.search(
                            rf'{re.escape(my_file)}\([+-]:(\d+)-(\d+)\|([0-9.]+)%\)',
                            other_details
                        )
                        
                        if other_match:
                            other_coverage_pct = float(other_match.group(3))
                            
                            # If I cover <=low_threshold% of the other range, remove me
                            if other_coverage_pct <= low_threshold:
                                return (True, other_file, my_coverage_pct, other_coverage_pct)
    
    return (False, None, 0.0, 0.0)

def resolve_by_priority(row, detailed_df):
    """
    Priority resolution for same-strand overlaps:
    1. More annotation tools (more reliable)
    2. If tie in tools, longer range (more complete)
    
    Returns: (should_remove, other_file, reason, my_value, other_value)
    """
    if row['n_overlapping_ranges'] == 0:
        return (False, None, None, 0, 0)
    
    if row['overlap_details'] == 'none':
        return (False, None, None, 0, 0)
    
    my_strand = row['strand']
    my_file = row['file']
    my_tools = row['n_annotation_tools']
    my_length = row['length']
    details = row['overlap_details'].split(';')
    
    for detail in details:
        # Parse: filename(strand:start-end|pct%)
        match = re.match(r'(.+?)\(([+-]):(\d+)-(\d+)\|([0-9.]+)%\)', detail.strip())
        if match:
            other_file = match.group(1)
            other_strand = match.group(2)
            
            # Only check same strand overlaps
            if other_strand != my_strand:
                continue
            
            other_row = detailed_df[detailed_df['file'] == other_file]
            if not other_row.empty:
                other_row = other_row.iloc[0]
                other_tools = other_row['n_annotation_tools']
                other_length = other_row['length']
                
                # Priority 1: More tools wins
                if other_tools > my_tools:
                    return (True, other_file, 'fewer_tools', my_tools, other_tools)
                elif other_tools == my_tools:
                    # Priority 2: Longer range wins
                    if other_length > my_length:
                        return (True, other_file, 'shorter_length', my_length, other_length)
    
    return (False, None, None, 0, 0)

def parse_args():
    parser = argparse.ArgumentParser(
        description='Filter singles based on multiple criteria'
    )
    parser.add_argument('--singles', required=True, help='singles.tsv file')
    parser.add_argument('--detailed', required=True, help='singles_detailed_analysis.tsv file')
    parser.add_argument('--groups', required=True, help='groups.tsv file')
    parser.add_argument('--coverage-threshold', type=float, default=90.0, 
                       help='Coverage threshold for group filtering (default: 90.0)')
    parser.add_argument('--overlap-threshold', type=float, default=30.0,
                       help='Overlap threshold for same-strand filtering (default: 30.0)')
    parser.add_argument('--asymmetric-high', type=float, default=80.0,
                       help='High threshold for asymmetric overlap detection (default: 80.0)')
    parser.add_argument('--asymmetric-low', type=float, default=45.0,
                       help='Low threshold for asymmetric overlap detection (default: 45.0)')
    parser.add_argument('--outdir', default='.', help='Output directory')
    return parser.parse_args()

def main():
    args = parse_args()
    
    print("Loading input files...")
    print(f"Coverage threshold: {args.coverage_threshold}%")
    print(f"Same-strand overlap threshold: {args.overlap_threshold}%")
    print(f"Asymmetric overlap thresholds: >{args.asymmetric_high}% covered, <={args.asymmetric_low}% covering")
    
    # Load files
    singles_df = pd.read_csv(args.singles, sep='\t')
    detailed_df = pd.read_csv(args.detailed, sep='\t')
    
    # Parse group ranges
    print("Parsing group ranges...")
    group_ranges = parse_group_ranges(args.groups)
    print(f"Found {len(group_ranges)} ranges in groups")
    
    # Get 5-tool ranges
    print("Identifying 5-tool ranges...")
    five_tool_ranges = get_five_tool_ranges(args.detailed)
    print(f"Found {len(five_tool_ranges)} ranges predicted by 5 tools")
    
    # Categories
    keep_five_tools = []
    keep_overlap_with_five = []
    keep_different_strand = []
    keep_same_strand_low_overlap = []
    keep_no_overlaps = []
    keep_by_priority = []  # NUEVA CATEGORÍA
    removed_covered = []
    removed_asymmetric = []
    removed_priority = []
    uncategorized = []
    
    # Track files to remove and winners
    files_to_remove = set()
    priority_winners = {}  # Track who won and why
    
    print("Filtering ranges...")
    
    # FIRST PASS: Identify ranges to remove (groups, asymmetric, and priority)
    for _, row in detailed_df.iterrows():
        filename = row['file']
        start = row['start']
        end = row['end']
        strand = row['strand']
        n_tools = row['n_annotation_tools']
        length = row['length']
        
        # Check if covered by groups (>= threshold%)
        is_covered, covering_file, coverage_pct = is_covered_by_groups(
            start, end, group_ranges, args.coverage_threshold
        )
        
        if is_covered:
            removed_covered.append({
                'file': filename,
                'start': start,
                'end': end,
                'strand': strand,
                'n_tools': n_tools,
                'length': length,
                'covered_by': covering_file,
                'coverage_percentage': round(coverage_pct, 2),
                'reason': f'covered_{args.coverage_threshold}%_by_group'
            })
            files_to_remove.add(filename)
            continue
        
        # Check for asymmetric overlap
        should_remove, other_file, my_cov, other_cov = check_asymmetric_overlap(
            row, detailed_df, args.asymmetric_high, args.asymmetric_low
        )
        
        if should_remove:
            removed_asymmetric.append({
                'file': filename,
                'start': start,
                'end': end,
                'strand': strand,
                'n_tools': n_tools,
                'length': length,
                'removed_in_favor_of': other_file,
                'my_coverage': round(my_cov, 2),
                'other_coverage': round(other_cov, 2),
                'reason': f'asymmetric_overlap_>{args.asymmetric_high}%_covered_<={args.asymmetric_low}%_covering'
            })
            files_to_remove.add(filename)
            continue
        
        # Check priority resolution (fewer tools or shorter length)
        should_remove, other_file, reason, my_val, other_val = resolve_by_priority(
            row, detailed_df
        )
        
        if should_remove:
            removed_priority.append({
                'file': filename,
                'start': start,
                'end': end,
                'strand': strand,
                'n_tools': n_tools,
                'length': length,
                'removed_in_favor_of': other_file,
                'reason': reason,
                'my_value': my_val,
                'other_value': other_val
            })
            files_to_remove.add(filename)
        else:
            # If not removed, check if it won a priority battle
            if reason is None and row['n_overlapping_ranges'] > 0:
                # Check if it has same-strand overlaps
                if not has_different_strand_overlaps(row) and not has_same_strand_low_overlap(row, args.overlap_threshold):
                    # This is a priority winner
                    priority_winners[filename] = {
                        'won_against': other_file if other_file else 'unknown',
                        'reason': 'priority_winner'
                    }
    
    print(f"Files marked for removal: {len(files_to_remove)}")
    print(f"Priority winners identified: {len(priority_winners)}")
    
    # SECOND PASS: Categorize remaining ranges
    for _, row in detailed_df.iterrows():
        filename = row['file']
        
        # Skip if already marked for removal
        if filename in files_to_remove:
            continue
        
        start = row['start']
        end = row['end']
        strand = row['strand']
        n_tools = row['n_annotation_tools']
        
        # Rule 1: Keep if 5 tools
        if n_tools == 5:
            keep_five_tools.append(row.to_dict())
            continue
        
        # Rule 2: Keep if overlaps ONLY with 5-tool ranges
        if overlaps_only_with_five_tools(row, five_tool_ranges):
            keep_overlap_with_five.append(row.to_dict())
            continue
        
        # Rule 3: Keep if has different strand overlaps
        if has_different_strand_overlaps(row):
            keep_different_strand.append(row.to_dict())
            continue
        
        # Rule 4: Keep if has same strand overlaps but ALL < threshold%
        if has_same_strand_low_overlap(row, args.overlap_threshold):
            keep_same_strand_low_overlap.append(row.to_dict())
            continue
        
        # Rule 5: Keep if has NO overlaps at all
        if has_no_overlaps(row):
            keep_no_overlaps.append(row.to_dict())
            continue
        
        # Rule 6: Keep if won priority battle (NEW)
        if filename in priority_winners:
            row_dict = row.to_dict()
            row_dict['priority_reason'] = priority_winners[filename]['reason']
            keep_by_priority.append(row_dict)
            continue
        
        # Otherwise, uncategorized
        uncategorized.append(row.to_dict())
    
    # Create output DataFrames
    print("\nCreating output files...")
    
    # 1. Ranges to keep (5 tools)
    if keep_five_tools:
        df_five = pd.DataFrame(keep_five_tools)
        df_five.to_csv(f"{args.outdir}/keep_five_tools.tsv", sep='\t', index=False)
        print(f"Ranges with 5 tools: {len(keep_five_tools)}")
    
    # 2. Ranges to keep (overlap with 5-tool ranges)
    if keep_overlap_with_five:
        df_overlap_five = pd.DataFrame(keep_overlap_with_five)
        df_overlap_five.to_csv(f"{args.outdir}/keep_overlap_with_five_tools.tsv", sep='\t', index=False)
        print(f"Ranges overlapping only with 5-tool ranges: {len(keep_overlap_with_five)}")
    
    # 3. Ranges to keep (different strand overlaps)
    if keep_different_strand:
        df_diff_strand = pd.DataFrame(keep_different_strand)
        df_diff_strand.to_csv(f"{args.outdir}/keep_different_strand_overlaps.tsv", sep='\t', index=False)
        print(f"Ranges with different strand overlaps: {len(keep_different_strand)}")
    
    # 4. Ranges to keep (same strand, low overlap)
    if keep_same_strand_low_overlap:
        df_same_low = pd.DataFrame(keep_same_strand_low_overlap)
        df_same_low.to_csv(f"{args.outdir}/keep_same_strand_low_overlap.tsv", sep='\t', index=False)
        print(f"Ranges with same strand, <{args.overlap_threshold}% overlap: {len(keep_same_strand_low_overlap)}")
    
    # 5. Ranges to keep (no overlaps)
    if keep_no_overlaps:
        df_no_overlaps = pd.DataFrame(keep_no_overlaps)
        df_no_overlaps.to_csv(f"{args.outdir}/keep_no_overlaps.tsv", sep='\t', index=False)
        print(f"Ranges with NO overlaps: {len(keep_no_overlaps)}")
    
    # 6. Ranges kept by priority (NEW)
    if keep_by_priority:
        df_priority = pd.DataFrame(keep_by_priority)
        df_priority.to_csv(f"{args.outdir}/keep_by_priority_resolution.tsv", sep='\t', index=False)
        print(f"Ranges kept by priority resolution: {len(keep_by_priority)}")
    
    # 7. All ranges to keep combined
    all_keep = (keep_five_tools + keep_overlap_with_five + keep_different_strand + 
                keep_same_strand_low_overlap + keep_no_overlaps + keep_by_priority)
    if all_keep:
        df_all_keep = pd.DataFrame(all_keep)
        df_all_keep.to_csv(f"{args.outdir}/all_ranges_to_keep.tsv", sep='\t', index=False)
        print(f"Total ranges to KEEP: {len(all_keep)}")
    
    # 8. Removed (covered by groups)
    if removed_covered:
        df_removed = pd.DataFrame(removed_covered)
        df_removed.to_csv(f"{args.outdir}/removed_covered_by_groups.tsv", sep='\t', index=False)
        print(f"Ranges REMOVED (covered >={args.coverage_threshold}% by groups): {len(removed_covered)}")
    
    # 9. Removed (asymmetric overlap)
    if removed_asymmetric:
        df_removed_asym = pd.DataFrame(removed_asymmetric)
        df_removed_asym.to_csv(f"{args.outdir}/removed_asymmetric_overlap.tsv", sep='\t', index=False)
        print(f"Ranges REMOVED (asymmetric overlap): {len(removed_asymmetric)}")
    
    # 10. Removed (priority: fewer tools or shorter)
    if removed_priority:
        df_removed_priority = pd.DataFrame(removed_priority)
        df_removed_priority.to_csv(f"{args.outdir}/removed_by_priority.tsv", sep='\t', index=False)
        print(f"Ranges REMOVED (priority: fewer tools or shorter): {len(removed_priority)}")
    
    # 11. Uncategorized
    if uncategorized:
        df_uncat = pd.DataFrame(uncategorized)
        df_uncat.to_csv(f"{args.outdir}/uncategorized_ranges.tsv", sep='\t', index=False)
        print(f"Uncategorized ranges: {len(uncategorized)}")
    
    # Summary
    print("\n=== FILTERING SUMMARY ===")
    print(f"Total input ranges: {len(detailed_df)}")
    print(f"Ranges with 5 tools: {len(keep_five_tools)}")
    print(f"Ranges overlapping only with 5-tool ranges: {len(keep_overlap_with_five)}")
    print(f"Ranges with different strand overlaps: {len(keep_different_strand)}")
    print(f"Ranges with same strand, <{args.overlap_threshold}% overlap: {len(keep_same_strand_low_overlap)}")
    print(f"Ranges with NO overlaps: {len(keep_no_overlaps)}")
    print(f"Ranges kept by priority resolution: {len(keep_by_priority)}")
    print(f"Total to KEEP: {len(all_keep)}")
    print(f"Removed (covered >={args.coverage_threshold}% by groups): {len(removed_covered)}")
    print(f"Removed (asymmetric overlap): {len(removed_asymmetric)}")
    print(f"Removed (priority resolution): {len(removed_priority)}")
    print(f"Uncategorized: {len(uncategorized)}")
    
    # Create summary file
    with open(f"{args.outdir}/filtering_summary.txt", 'w') as f:
        f.write("=== SINGLES FILTERING SUMMARY ===\n\n")
        f.write(f"Coverage threshold: {args.coverage_threshold}%\n")
        f.write(f"Same-strand overlap threshold: {args.overlap_threshold}%\n")
        f.write(f"Asymmetric overlap thresholds: >{args.asymmetric_high}% covered, <={args.asymmetric_low}% covering\n\n")
        f.write(f"Total input ranges: {len(detailed_df)}\n\n")
        f.write("RANGES TO KEEP:\n")
        f.write(f"  - Predicted by 5 tools: {len(keep_five_tools)}\n")
        f.write(f"  - Overlapping only with 5-tool ranges: {len(keep_overlap_with_five)}\n")
        f.write(f"  - With different strand overlaps: {len(keep_different_strand)}\n")
        f.write(f"  - With same strand, <{args.overlap_threshold}% overlap: {len(keep_same_strand_low_overlap)}\n")
        f.write(f"  - With NO overlaps at all: {len(keep_no_overlaps)}\n")
        f.write(f"  - Kept by priority resolution: {len(keep_by_priority)}\n")
        f.write(f"  TOTAL TO KEEP: {len(all_keep)}\n\n")
        f.write("RANGES REMOVED:\n")
        f.write(f"  - Covered >={args.coverage_threshold}% by groups: {len(removed_covered)}\n")
        f.write(f"  - Asymmetric overlap (>{args.asymmetric_high}% covered, <={args.asymmetric_low}% covering): {len(removed_asymmetric)}\n")
        f.write(f"  - Priority resolution (fewer tools or shorter): {len(removed_priority)}\n")
        f.write(f"  TOTAL REMOVED: {len(removed_covered) + len(removed_asymmetric) + len(removed_priority)}\n\n")
        f.write(f"UNCATEGORIZED: {len(uncategorized)}\n")

if __name__ == "__main__":
    main()
