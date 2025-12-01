#!/usr/bin/env python3

"""
Update merged_final_ranges files with identified frames.
"""

import sys
import argparse
import pandas as pd
import re

def load_identified_frames(frames_file):
    """
    Load identified frames into a dictionary.
    Creates keys for both original and +1 start positions to handle bedtools offset.
    """
    df = pd.read_csv(frames_file, sep='\t')
    frames_dict = {}
    
    for _, row in df.iterrows():
        # Parse range: start-end(strand)
        match = re.match(r'(\d+)-(\d+)\(([+-])\)', row['range'])
        if match:
            start = int(match.group(1))
            end = int(match.group(2))
            
            # Create entries for both possible start positions
            # Original coordinates (from BED name field)
            key1 = (start, end)
            # Adjusted coordinates (+1 start, common after bedtools)
            key2 = (start + 1, end)
            
            frame_info = {
                'frame': int(row['frame']),
                'status': row['status']
            }
            
            frames_dict[key1] = frame_info
            frames_dict[key2] = frame_info
    
    return frames_dict

def update_tsv(input_tsv, frames_dict, output_tsv):
    """Update TSV file with identified frames."""
    df = pd.read_csv(input_tsv, sep='\t')
    
    updated_count = 0
    not_found = []
    
    for idx, row in df.iterrows():
        if row['frame'] == 'unknown':
            key = (int(row['start']), int(row['end']))
            if key in frames_dict:
                df.at[idx, 'frame'] = frames_dict[key]['frame']
                updated_count += 1
            else:
                not_found.append(f"{row['start']}-{row['end']}")
    
    df.to_csv(output_tsv, sep='\t', index=False)
    
    if not_found:
        print(f"Warning: {len(not_found)} ranges not found in frames dict:", file=sys.stderr)
        for nf in not_found[:5]:
            print(f"  {nf}", file=sys.stderr)
        if len(not_found) > 5:
            print(f"  ... and {len(not_found)-5} more", file=sys.stderr)
    
    return updated_count

def update_bed(input_bed, frames_dict, output_bed, chrom_name):
    """Update BED file with identified frames."""
    updated_count = 0
    not_found = []
    
    with open(input_bed, 'r') as infile, open(output_bed, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) < 6:
                outfile.write(line)
                continue
            
            chrom, start_0based, end_pos, name, score, strand = fields
            
            # Parse name: start-end_strand_frame
            match = re.match(r'(\d+)-(\d+)_([+.-])_(unk|\d+)', name)
            if match and match.group(4) == 'unk':
                start = int(match.group(1))
                end = int(match.group(2))
                
                # Try both possible keys
                key1 = (start, end)
                key2 = (start + 1, end)
                
                frame_info = None
                if key1 in frames_dict:
                    frame_info = frames_dict[key1]
                elif key2 in frames_dict:
                    frame_info = frames_dict[key2]
                
                if frame_info:
                    new_frame = frame_info['frame']
                    new_name = f"{start}-{end}_{strand}_{new_frame}"
                    fields[3] = new_name
                    updated_count += 1
                else:
                    not_found.append(f"{start}-{end}")
            
            outfile.write('\t'.join(fields) + '\n')
    
    if not_found:
        print(f"Warning: {len(not_found)} ranges not found in frames dict (BED):", file=sys.stderr)
        for nf in not_found[:5]:
            print(f"  {nf}", file=sys.stderr)
        if len(not_found) > 5:
            print(f"  ... and {len(not_found)-5} more", file=sys.stderr)
    
    return updated_count

def parse_args():
    parser = argparse.ArgumentParser(
        description='Update merged ranges files with identified frames'
    )
    parser.add_argument('--frames', required=True,
                       help='identified_frames.tsv file')
    parser.add_argument('--input-tsv', required=True,
                       help='Input merged_final_ranges.tsv')
    parser.add_argument('--input-bed', required=True,
                       help='Input merged_final_ranges.bed')
    parser.add_argument('--output-tsv', required=True,
                       help='Output updated TSV')
    parser.add_argument('--output-bed', required=True,
                       help='Output updated BED')
    parser.add_argument('--chrom-name', default='Lokiarchaeum_B36',
                       help='Chromosome name for BED file')
    return parser.parse_args()

def main():
    args = parse_args()
    
    print("Loading identified frames...")
    frames_dict = load_identified_frames(args.frames)
    # Divide by 2 because we create duplicate entries
    unique_frames = len(frames_dict) // 2
    print(f"  Loaded {unique_frames} frame assignments (with coordinate variants)")
    
    print("\nUpdating TSV file...")
    tsv_updated = update_tsv(args.input_tsv, frames_dict, args.output_tsv)
    print(f"  Updated {tsv_updated} ranges in TSV")
    
    print("\nUpdating BED file...")
    bed_updated = update_bed(args.input_bed, frames_dict, args.output_bed, args.chrom_name)
    print(f"  Updated {bed_updated} ranges in BED")
    
    # Summary statistics
    print("\n=== UPDATE SUMMARY ===")
    print(f"Total frame assignments available: {unique_frames}")
    print(f"Ranges updated in TSV: {tsv_updated}")
    print(f"Ranges updated in BED: {bed_updated}")
    
    # Verify all unknowns were updated
    df_updated = pd.read_csv(args.output_tsv, sep='\t')
    remaining_unknown = (df_updated['frame'] == 'unknown').sum()
    
    if remaining_unknown > 0:
        print(f"\nWARNING: {remaining_unknown} ranges still have unknown frame!")
    else:
        print(f"\n✓ All frames successfully identified and updated!")
    
    print(f"\nUpdated files saved:")
    print(f"  TSV: {args.output_tsv}")
    print(f"  BED: {args.output_bed}")

if __name__ == "__main__":
    main()

