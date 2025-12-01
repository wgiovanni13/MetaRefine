#!/usr/bin/env python3
"""
Create a simple summary of frameshift analysis results
"""

import sys
from collections import Counter

def main():
    if len(sys.argv) < 3:
        print("Usage: summarize_frameshifts.py <input.tsv> <output.txt>", file=sys.stderr)
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Read frameshifts
    frameshifts = []
    with open(input_file, 'r') as f:
        header = f.readline().strip().split('\t')
        
        # Debug: print header
        print(f"[INFO] Reading {input_file}", file=sys.stderr)
        print(f"[INFO] Header has {len(header)} columns", file=sys.stderr)
        
        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line:
                continue
                
            fields = line.split('\t')
            
            if len(fields) < 26:
                print(f"[WARNING] Line {line_num} has only {len(fields)} columns, skipping", file=sys.stderr)
                continue
            
            try:
                frameshifts.append({
                    'gene1': fields[0],
                    'gene2': fields[1],
                    'ref_id': fields[7],
                    'target_desc': fields[8],
                    'avg_identity': float(fields[12]),
                    'total_cov': float(fields[21]),
                    'confidence': fields[25]
                })
            except (ValueError, IndexError) as e:
                print(f"[WARNING] Error parsing line {line_num}: {e}", file=sys.stderr)
                continue
    
    print(f"[INFO] Successfully parsed {len(frameshifts)} frameshifts", file=sys.stderr)
    
    # Generate summary
    with open(output_file, 'w') as out:
        out.write("=" * 80 + "\n")
        out.write("FRAMESHIFT ANALYSIS SUMMARY\n")
        out.write("=" * 80 + "\n\n")
        
        out.write(f"Total frameshift candidates detected: {len(frameshifts)}\n\n")
        
        if len(frameshifts) == 0:
            out.write("No frameshifts detected.\n")
            out.write("\nNote: Only characterized proteins are reported.\n")
            out.write("Uncharacterized/hypothetical proteins are filtered out.\n")
            return
        
        # Confidence breakdown
        confidence_counts = Counter(fs['confidence'] for fs in frameshifts)
        out.write("Confidence levels:\n")
        for conf in ['HIGH', 'MEDIUM', 'LOW']:
            count = confidence_counts.get(conf, 0)
            pct = (count / len(frameshifts) * 100) if len(frameshifts) > 0 else 0
            out.write(f"  {conf}: {count} ({pct:.1f}%)\n")
        out.write("\n")
        
        # Top proteins
        protein_counts = Counter(fs['ref_id'] for fs in frameshifts)
        out.write("Top 10 proteins with most frameshifts:\n")
        for protein, count in protein_counts.most_common(10):
            # Get description from first occurrence
            desc = next(fs['target_desc'] for fs in frameshifts if fs['ref_id'] == protein)
            out.write(f"  {count}x - {protein}\n")
            out.write(f"       {desc[:70]}...\n" if len(desc) > 70 else f"       {desc}\n")
        out.write("\n")
        
        # High-confidence frameshifts
        high_conf = [fs for fs in frameshifts if fs['confidence'] == 'HIGH']
        out.write(f"High-confidence frameshifts: {len(high_conf)}\n")
        out.write("-" * 80 + "\n\n")
        
        if high_conf:
            # Sort by avg_identity
            high_conf.sort(key=lambda x: x['avg_identity'], reverse=True)
            
            out.write("Top 10 high-confidence frameshifts (by identity):\n\n")
            for i, fs in enumerate(high_conf[:10], 1):
                out.write(f"{i}. Genes: {fs['gene1']} + {fs['gene2']}\n")
                out.write(f"   Reference: {fs['ref_id']}\n")
                desc_short = fs['target_desc'][:65] + "..." if len(fs['target_desc']) > 65 else fs['target_desc']
                out.write(f"   Description: {desc_short}\n")
                out.write(f"   Identity: {fs['avg_identity']:.1f}%  |  Coverage: {fs['total_cov']:.1f}%\n\n")
        
        out.write("\n" + "=" * 80 + "\n")
        out.write(f"Total: {len(frameshifts)} characterized frameshifts detected\n")
        out.write("(Uncharacterized/hypothetical proteins excluded)\n")
        out.write("Full details in: all_frameshifts.tsv\n")
        out.write("=" * 80 + "\n")
    
    print(f"[INFO] Summary written to {output_file}", file=sys.stderr)

if __name__ == '__main__':
    main()

