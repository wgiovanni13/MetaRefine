#!/usr/bin/env python3

"""
Identify correct reading frames from 6-frame translations.
Correct frame = exactly 1 stop codon at the end of the sequence.
"""

import sys
import argparse
import re
from collections import defaultdict

def parse_fasta(fasta_file):
    """Parse FASTA file and return dict of sequences."""
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_header:
            sequences[current_header] = ''.join(current_seq)
    
    return sequences

def count_stops(sequence):
    """
    Count internal stop codons (excluding terminal).
    Returns: (internal_stops, has_terminal_stop)
    """
    has_terminal = sequence.endswith('*')
    
    if has_terminal:
        body = sequence[:-1]  # Remove terminal stop
    else:
        body = sequence
    
    internal_stops = body.count('*')
    
    return internal_stops, has_terminal

def parse_header(header):
    """
    Parse header to extract range, strand, and frame.
    Format: start-end(strand)_frame
    Example: 687-1620(+)_1
    Returns: (range_key, strand, frame)
    """
    # Pattern: digits-digits(+/-)_digits
    match = re.match(r'(\d+)-(\d+)\(([+-])\)_(\d+)', header)
    if match:
        start = match.group(1)
        end = match.group(2)
        strand = match.group(3)
        frame = match.group(4)
        range_key = f"{start}-{end}({strand})"
        return range_key, strand, int(frame)
    return None, None, None

def identify_correct_frames(sequences):
    """
    Identify correct frame for each range.
    Correct frame = exactly 0 internal stops + 1 terminal stop.
    """
    # Group sequences by range
    range_frames = defaultdict(list)
    
    for header, seq in sequences.items():
        range_key, strand, frame = parse_header(header)
        if range_key is None:
            print(f"Warning: Could not parse header: {header}", file=sys.stderr)
            continue
        
        internal_stops, has_terminal = count_stops(seq)
        
        range_frames[range_key].append({
            'frame': frame,
            'strand': strand,
            'internal_stops': internal_stops,
            'has_terminal': has_terminal,
            'seq': seq,
            'header': header
        })
    
    # Identify correct frame for each range
    correct_frames = {}
    ambiguous_ranges = []
    no_valid_frame = []
    
    for range_key, frames in range_frames.items():
        # Find frames with exactly 0 internal stops and 1 terminal stop
        valid_frames = [f for f in frames 
                       if f['internal_stops'] == 0 and f['has_terminal']]
        
        if len(valid_frames) == 1:
            # Perfect: exactly one valid frame
            correct_frames[range_key] = valid_frames[0]
        elif len(valid_frames) > 1:
            # Ambiguous: multiple valid frames
            ambiguous_ranges.append((range_key, valid_frames))
            # Choose the first one as default
            correct_frames[range_key] = valid_frames[0]
        else:
            # No valid frame found
            no_valid_frame.append(range_key)
            # Choose frame with minimum internal stops
            best = min(frames, key=lambda x: x['internal_stops'])
            correct_frames[range_key] = best
    
    return correct_frames, ambiguous_ranges, no_valid_frame

def parse_args():
    parser = argparse.ArgumentParser(
        description='Identify correct reading frames from 6-frame translations'
    )
    parser.add_argument('--faa', required=True,
                       help='6-frame translated FASTA file')
    parser.add_argument('--output', default='identified_frames.tsv',
                       help='Output TSV with identified frames')
    return parser.parse_args()

def main():
    args = parse_args()
    
    print("Loading 6-frame translations...")
    sequences = parse_fasta(args.faa)
    print(f"  Loaded {len(sequences)} sequences")
    
    print("\nIdentifying correct frames...")
    correct_frames, ambiguous, no_valid = identify_correct_frames(sequences)
    
    print(f"  Ranges with unique correct frame: {len(correct_frames) - len(ambiguous) - len(no_valid)}")
    print(f"  Ranges with ambiguous frames: {len(ambiguous)}")
    print(f"  Ranges with no valid frame: {len(no_valid)}")
    
    # Write output
    with open(args.output, 'w') as out:
        out.write("range\tstrand\tframe\tinternal_stops\thas_terminal_stop\tstatus\n")
        
        for range_key, info in sorted(correct_frames.items()):
            # Determine status
            if range_key in [r for r, _ in ambiguous]:
                status = "ambiguous"
            elif range_key in no_valid:
                status = "no_valid_frame"
            else:
                status = "unique"
            
            out.write(f"{range_key}\t{info['strand']}\t{info['frame']}\t"
                     f"{info['internal_stops']}\t{info['has_terminal']}\t{status}\n")
    
    print(f"\nIdentified frames saved to: {args.output}")
    
    # Print warnings for ambiguous and no valid frames
    if ambiguous:
        print(f"\nWARNING: {len(ambiguous)} ranges have multiple valid frames:")
        for range_key, frames in ambiguous[:5]:  # Show first 5
            frame_list = [str(f['frame']) for f in frames]
            print(f"  {range_key}: frames {', '.join(frame_list)}")
        if len(ambiguous) > 5:
            print(f"  ... and {len(ambiguous) - 5} more")
    
    if no_valid:
        print(f"\nWARNING: {len(no_valid)} ranges have no valid frame (using best guess):")
        for range_key in no_valid[:5]:  # Show first 5
            info = correct_frames[range_key]
            print(f"  {range_key}: frame {info['frame']} (internal stops: {info['internal_stops']})")
        if len(no_valid) > 5:
            print(f"  ... and {len(no_valid) - 5} more")

if __name__ == "__main__":
    main()

