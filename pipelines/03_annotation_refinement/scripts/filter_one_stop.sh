#!/usr/bin/env bash

#SBATCH --job-name=one_stop_codon
#SBATCH --partition=c
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=one_stop_%j.out
#SBATCH --error=one_stop_%j.err
#SBATCH --mail-user=wagner.guzman@gmi.oeaw.ac.at
#SBATCH --mail-type=FAIL

# Filter sequences to keep only those with:
# 1. Exactly ONE frame (out of 6) has a terminal stop codon
# 2. That frame has NO internal stop codons (only 1 stop total, at the end)

FASTA="$1"   # cat3_ranges_extracted.6frames.archaea.faa
LIST="$2"    # ranges_with_exactly1_terminal_stop.txt
OUT="$3"     # cat3_ranges_extracted.6frames.archaea.oneStopStripped.faa

echo "=== Filtering for perfect ORFs (1 terminal stop, no internal stops) ==="
echo "Input: $FASTA"
echo "Output list: $LIST"
echo "Output FASTA: $OUT"
echo ""

# Step 1: Find sequences with exactly 1 stop at the end and no internal stops

awk '
BEGIN{ 
    id=""; seq="" 
    frames_with_terminal = 0
    frames_with_only_terminal = 0
}

function process(){
    if(id=="") return
    gsub(/[ \t\r\n]/,"",seq)
    if(seq=="") return
    
    # Extract base ID (remove _1, _2, etc.)
    base = id
    sub(/_[1-6]$/,"",base)
    
    # Count stop codons
    total_stops = gsub(/\*/,"*",seq)  # Total * in sequence
    has_terminal = (substr(seq, length(seq), 1) == "*")
    
    if(has_terminal) {
        frames_terminal[base]++
        
        # Check if this is the ONLY stop (no internal stops)
        if(total_stops == 1) {
            frames_perfect[base]++
            frame_id[base] = id  # Store the frame ID
        }
    }
    
    seq = ""
}

/^>/{
    process()
    id = $0
    sub(/^>/,"",id)
    split(id, a, /[ \t]/)
    id = a[1]
    next
}

{ 
    gsub(/[ \t\r]/,"")
    seq = seq $0 
}

END{
    process()
    
    # Select base IDs where:
    # - Exactly 1 frame has terminal stop
    # - That frame has no internal stops
    for (b in frames_terminal) {
        if (frames_terminal[b] == 1 && frames_perfect[b] == 1) {
            print frame_id[b]
            count++
        }
    }
    
    printf("Total base sequences: %d\n", length(frames_terminal)) > "/dev/stderr"
    printf("Perfect ORFs (1 terminal stop, no internal): %d\n", count) > "/dev/stderr"
}
' "$FASTA" | sort -u > "$LIST"

echo "Perfect ORFs found: $(wc -l < "$LIST")"
echo ""

# Step 2: Extract those sequences and strip the final stop codon

awk -v LIST="$LIST" '
BEGIN{
    # Load list of IDs to keep
    while((getline l < LIST) > 0) { 
        gsub(/\r/,"",l)
        keep[l] = 1
    }
    close(LIST)
    id = ""
    seq = ""
    kept = 0
    skipped = 0
}

function emit(){
    if(id == "") return
    gsub(/[ \t\r\n]/,"", seq)
    if(seq == "") return
    
    if(id in keep) {
        # Verify it has exactly 1 stop at the end
        total_stops = gsub(/\*/,"*",seq)
        has_terminal = (substr(seq, length(seq), 1) == "*")
        
        if(has_terminal && total_stops == 1) {
            # Strip terminal stop
            seq = substr(seq, 1, length(seq) - 1)
            
            # Write output
            print ">" id
            for(i = 1; i <= length(seq); i += 60) {
                print substr(seq, i, 60)
            }
            kept++
        } else {
            printf("WARNING: %s has %d stops but was in list\n", id, total_stops) > "/dev/stderr"
            skipped++
        }
    }
    seq = ""
}

/^>/{
    emit()
    id = $0
    sub(/^>/,"", id)
    split(id, a, /[ \t]/)
    id = a[1]
    next
}

{ 
    gsub(/[ \t\r]/,"")
    seq = seq $0 
}

END{
    emit()
    printf("Sequences kept: %d\n", kept) > "/dev/stderr"
    printf("Sequences skipped (quality issues): %d\n", skipped) > "/dev/stderr"
}
' "$FASTA" > "$OUT"

echo ""
echo "=== SUMMARY ==="
echo "List file: $$LIST ($$ (wc -l < "$LIST") IDs)"
echo "Output FASTA: $$OUT ($$ (grep -c '^>' "$OUT") sequences)"
echo ""

# Verification: check for any remaining internal stops

internal_stops=$(awk '/^>/{next} {gsub(/\*/,"*"); if(gsub(/\*/,"",$0)>0) count++} END{print count+0}' "$OUT")
if [ "$internal_stops" -gt 0 ]; then
    echo "WARNING: Found $internal_stops sequences with internal stops!"
else
    echo "✓ All sequences are clean (no internal stops)"
fi
