#!/usr/bin/env python3
"""
Extract the correct frame from 6-frame translated sequences.
Final quality control and cleanup step.

Strategy: 
1. Select the translation with exactly ONE stop codon at the END (perfect ORF)
2. REJECT any sequences with internal stop codons
3. Remove terminal stop codon from all selected sequences

IMPORTANT: Input nucleotide sequences were extracted with bedtools getfasta -s,
which already applied reverse complement for minus strand sequences.
"""

import argparse
import re
import sys
from collections import defaultdict

def read_fasta(filename):
    """Read FASTA file and return dict of id -> sequence."""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def count_stop_codons(sequence):
    """
    Count stop codons in protein sequence.
    Stop codons are represented as: * (asterisk)
    
    Returns: (total_stops, internal_stops, has_terminal_stop)
    """
    total_stops = sequence.count('*')
    has_terminal_stop = sequence.endswith('*')
    
    if has_terminal_stop:
        internal_stops = total_stops - 1
    else:
        internal_stops = total_stops
    
    return total_stops, internal_stops, has_terminal_stop

def get_base_id(header):
    """
    Extract base ID without the transeq frame suffix.
    
    Input: 176-545(-)_1
    Output: 176-545(-)
    """
    base = re.sub(r'_[1-6]$', '', header)
    return base

def select_best_frame(frame_group):
    """
    Select the best frame from a group of 6 translations.
    
    STRICT CRITERIA (quality control):
    1. MUST have exactly 1 stop codon at the end
    2. MUST NOT have any internal stop codons
    3. If multiple frames pass, select longest
    
    Returns: best_data dict or None if no valid frame found
    """
    valid_candidates = []
    
    for data in frame_group:
        seq = data['seq']
        total_stops, internal_stops, has_terminal = count_stop_codons(seq)
        
        # STRICT: Only accept perfect ORFs (1 terminal stop, 0 internal stops)
        if total_stops == 1 and has_terminal and internal_stops == 0:
            valid_candidates.append({
                'data': data,
                'length': len(seq),
                'total_stops': total_stops,
                'internal_stops': internal_stops,
                'has_terminal': has_terminal
            })
    
    if not valid_candidates:
        return None
    
    # Sort by length (longest first)
    valid_candidates.sort(key=lambda x: -x['length'])
    
    return valid_candidates[0]

def extract_correct_frames(transeq_fasta, output_fasta, stats_file, rejected_file):
    """
    Extract the correct frame for each sequence group.
    STRICT quality control: only perfect ORFs pass.
    """
    print(f"\n=== EXTRACTING CORRECT FRAMES (STRICT QC) ===")
    print(f"Input: {transeq_fasta}")
    print(f"Strategy: Select frame with exactly 1 terminal stop, NO internal stops")
    print()
    
    # Read all sequences
    try:
        all_seqs = read_fasta(transeq_fasta)
        print(f"Total translated sequences: {len(all_seqs):,}")
    except Exception as e:
        print(f"ERROR reading FASTA: {e}")
        sys.exit(1)
    
    # Group sequences by base ID
    grouped = defaultdict(list)
    
    for header, seq in all_seqs.items():
        # Extract transeq frame number
        match = re.search(r'_([1-6])$', header)
        if not match:
            print(f"WARNING: No frame suffix in header: {header}")
            continue
        
        frame_num = int(match.group(1))
        base_id = get_base_id(header)
        
        grouped[base_id].append({
            'header': header,
            'seq': seq,
            'frame': frame_num,
            'base_id': base_id
        })
    
    print(f"Unique sequences (groups): {len(grouped):,}")
    print()
    
    # Select best frame for each group
    selected = []
    rejected = []
    stats_data = []
    
    perfect_count = 0
    rejected_count = 0
    
    for base_id, frames in grouped.items():
        if len(frames) != 6:
            print(f"WARNING: {base_id} has {len(frames)} frames (expected 6)")
        
        best = select_best_frame(frames)
        
        if best:
            # Perfect ORF found
            selected.append(best['data'])
            perfect_count += 1
            
            stats_data.append({
                'base_id': base_id,
                'frame': best['data']['frame'],
                'length': best['length'],
                'total_stops': best['total_stops'],
                'internal_stops': best['internal_stops'],
                'has_terminal': best['has_terminal'],
                'status': 'accepted'
            })
        else:
            # No valid frame - all have internal stops or no terminal stop
            rejected_count += 1
            
            # Find the "best" among bad options for reporting
            all_stats = []
            for frame_data in frames:
                seq = frame_data['seq']
                total_stops, internal_stops, has_terminal = count_stop_codons(seq)
                all_stats.append({
                    'frame': frame_data['frame'],
                    'total_stops': total_stops,
                    'internal_stops': internal_stops,
                    'has_terminal': has_terminal
                })
            
            rejected.append({
                'base_id': base_id,
                'frames_stats': all_stats
            })
            
            # Add to stats as rejected
            stats_data.append({
                'base_id': base_id,
                'frame': 'NA',
                'length': 'NA',
                'total_stops': 'NA',
                'internal_stops': 'NA',
                'has_terminal': 'NA',
                'status': 'rejected_internal_stops'
            })
    
    print(f"=== QUALITY CONTROL SUMMARY ===")
    print(f"Perfect ORFs (accepted):      {perfect_count:>6,} ({perfect_count/(perfect_count+rejected_count)*100:>5.1f}%)")
    print(f"Rejected (internal stops):    {rejected_count:>6,} ({rejected_count/(perfect_count+rejected_count)*100:>5.1f}%)")
    print(f"Total:                        {perfect_count+rejected_count:>6,}")
    print()
    
    # Write rejected sequences report
    with open(rejected_file, 'w') as rej_fh:
        rej_fh.write("# Sequences rejected due to internal stop codons\n")
        rej_fh.write("# All 6 frames had either internal stops or no terminal stop\n\n")
        
        for rej in rejected:
            rej_fh.write(f"Sequence: {rej['base_id']}\n")
            for frame_stat in rej['frames_stats']:
                rej_fh.write(f"  Frame {frame_stat['frame']}: "
                           f"total_stops={frame_stat['total_stops']}, "
                           f"internal_stops={frame_stat['internal_stops']}, "
                           f"has_terminal={frame_stat['has_terminal']}\n")
            rej_fh.write("\n")
    
    # Write output FASTA (with terminal stops REMOVED)
    total = 0
    with open(output_fasta, 'w') as out_fh, open(stats_file, 'w') as stats_fh:
        # Write stats header
        stats_fh.write("base_id\tframe\tlength_with_stop\tlength_without_stop\t"
                      "total_stops\tinternal_stops\thas_terminal_stop\tstatus\n")
        
        # Write sequences
        for data in selected:
            total += 1
            
            seq = data['seq']
            
            # Remove terminal stop codon
            if seq.endswith('*'):
                seq_clean = seq[:-1]
            else:
                seq_clean = seq
                print(f"WARNING: {data['header']} doesn't end with * but was selected?")
            
            # Write FASTA header
            out_fh.write(f">{data['header']}\n")
            
            # Write sequence (60 chars per line)
            for i in range(0, len(seq_clean), 60):
                out_fh.write(f"{seq_clean[i:i+60]}\n")
        
        # Write stats
        for stat in stats_data:
            if stat['status'] == 'accepted':
                length_with = stat['length']
                length_without = length_with - 1  # removed terminal *
                stats_fh.write(
                    f"{stat['base_id']}\t{stat['frame']}\t{length_with}\t{length_without}\t"
                    f"{stat['total_stops']}\t{stat['internal_stops']}\t"
                    f"{stat['has_terminal']}\t{stat['status']}\n"
                )
            else:
                stats_fh.write(
                    f"{stat['base_id']}\t{stat['frame']}\t{stat['length']}\t"
                    f"{stat['length']}\t{stat['total_stops']}\t{stat['internal_stops']}\t"
                    f"{stat['has_terminal']}\t{stat['status']}\n"
                )
    
    print(f"=== OUTPUT ===")
    print(f"Proteins (clean, no stops): {output_fasta}")
    print(f"Stats: {stats_file}")
    print(f"Rejected report: {rejected_file}")
    print(f"Total sequences written: {total:,}")
    print()
    
    # Final verification
    print(f"=== FINAL VERIFICATION ===")
    with open(output_fasta) as f:
        clean_seqs = read_fasta(output_fasta)
    
    stops_found = 0
    for seq_id, seq in clean_seqs.items():
        if '*' in seq:
            stops_found += 1
            print(f"ERROR: Stop codon found in {seq_id}")
    
    if stops_found == 0:
        print(f"✓ All {len(clean_seqs):,} sequences are clean (no stop codons)")
    else:
        print(f"✗ ERROR: Found {stops_found} sequences with stop codons!")
        sys.exit(1)
    
    # Show examples
    print(f"\n=== EXAMPLES ===")
    print("First 3 sequences:")
    for i, (seq_id, seq) in enumerate(list(clean_seqs.items())[:3]):
        print(f"  {seq_id}: length {len(seq)} aa")
    
    if rejected_count > 0:
        print(f"\nRejected sequences (see {rejected_file} for details):")
        for rej in rejected[:3]:
            print(f"  {rej['base_id']}")
    
    if total == 0:
        print("\nERROR: No sequences passed quality control!")
        sys.exit(1)

def parse_args():
    parser = argparse.ArgumentParser(
        description='Extract correct frame with strict quality control',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
STRICT Selection criteria:
  - MUST have exactly 1 stop codon at the END
  - MUST NOT have any internal stop codons
  - Terminal stop codon is REMOVED from output
  - Sequences failing QC are rejected and reported

Input sequences were extracted with bedtools -s, so reverse
complement was already applied for minus strand sequences.
        """
    )
    parser.add_argument('--input', required=True,
                       help='Input FASTA (all 6 frames from transeq)')
    parser.add_argument('--output', default='final_proteins.faa',
                       help='Output FASTA (clean proteins, no stops)')
    parser.add_argument('--stats', default='frame_selection_stats.tsv',
                       help='Output statistics file')
    parser.add_argument('--rejected', default='rejected_sequences.txt',
                       help='Report of rejected sequences')
    return parser.parse_args()

def main():
    args = parse_args()
    
    try:
        extract_correct_frames(args.input, args.output, args.stats, args.rejected)
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

