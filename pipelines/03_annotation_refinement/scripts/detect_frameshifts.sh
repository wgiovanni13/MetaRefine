#!/bin/bash

INPUT="$1"
OUTPUT="$2"

echo -e "gene1\tgene2\tgene1_coords\tgene2_coords\tgene1_strand\tgene2_strand\tdistance_bp\tref_id\tref_length\tgene1_identity\tgene2_identity\tavg_identity\tgene1_evalue\tgene2_evalue\tgene1_ref_start\tgene1_ref_end\tgene2_ref_start\tgene2_ref_end\tgene1_cov\tgene2_cov\ttotal_cov\tgap_in_ref\toverlap_aa\tref_order\tconfidence" > "$OUTPUT"

zcat "$INPUT" | awk -F'\t' '
BEGIN {
    min_identity = 30
    max_evalue = 0.001
    max_overlap = 0.3
}

{
    gene = $1
    target = $2
    pident = $3
    evalue = $11
    t_start = $9
    t_end = $10
    t_len = $14
    
    if (pident >= min_identity && evalue <= max_evalue) {
        key = gene "|" target
        if (!(key in stored)) {
            genes[gene] = 1
            targets[gene][target]["pident"] = pident
            targets[gene][target]["evalue"] = evalue
            targets[gene][target]["t_start"] = t_start
            targets[gene][target]["t_end"] = t_end
            targets[gene][target]["t_len"] = t_len
            stored[key] = 1
        }
    }
}

END {
    n = 0
    for (gene in genes) {
        # Parse: 674504-675272(-)_1
        # Split by dash first
        split(gene, p1, "-")
        start = p1[1] + 1  # BED to genomic
        
        # Now parse the rest: 675272(-)_1
        rest = p1[2]
        
        # Find the opening paren position
        paren_pos = index(rest, "(")
        if (paren_pos > 0) {
            end = substr(rest, 1, paren_pos - 1)
            
            # Extract strand: character after (
            strand = substr(rest, paren_pos + 1, 1)
            
            # Extract frame: character after )_
            underscore_pos = index(rest, "_")
            if (underscore_pos > 0) {
                frame = substr(rest, underscore_pos + 1)
                
                sorted_genes[++n] = gene
                gene_start[gene] = start + 0  # Force numeric
                gene_end[gene] = end + 0
                gene_strand[gene] = strand
                gene_frame[gene] = frame
            }
        }
    }
    
    # Sort by position
    for (i = 1; i <= n; i++) {
        for (j = i+1; j <= n; j++) {
            if (gene_start[sorted_genes[i]] > gene_start[sorted_genes[j]]) {
                temp = sorted_genes[i]
                sorted_genes[i] = sorted_genes[j]
                sorted_genes[j] = temp
            }
        }
    }
    
    print "Analyzing", n, "valid genes..." > "/dev/stderr"
    
    frameshifts = 0
    for (i = 1; i < n; i++) {
        gene1 = sorted_genes[i]
        gene2 = sorted_genes[i+1]
        
        s1 = gene_start[gene1]
        e1 = gene_end[gene1]
        strand1 = gene_strand[gene1]
        
        s2 = gene_start[gene2]
        e2 = gene_end[gene2]
        strand2 = gene_strand[gene2]
        
        distance = s2 - e1 - 1
        
        for (t in targets[gene1]) {
            if (t in targets[gene2]) {
                pident1 = targets[gene1][t]["pident"]
                pident2 = targets[gene2][t]["pident"]
                evalue1 = targets[gene1][t]["evalue"]
                evalue2 = targets[gene2][t]["evalue"]
                t_start1 = targets[gene1][t]["t_start"]
                t_end1 = targets[gene1][t]["t_end"]
                t_start2 = targets[gene2][t]["t_start"]
                t_end2 = targets[gene2][t]["t_end"]
                t_len = targets[gene1][t]["t_len"]
                
                overlap_start = (t_start1 > t_start2) ? t_start1 : t_start2
                overlap_end = (t_end1 < t_end2) ? t_end1 : t_end2
                overlap = (overlap_end >= overlap_start) ? (overlap_end - overlap_start + 1) : 0
                
                len1 = t_end1 - t_start1 + 1
                len2 = t_end2 - t_start2 + 1
                
                ovlp_frac1 = (len1 > 0) ? overlap / len1 : 0
                ovlp_frac2 = (len2 > 0) ? overlap / len2 : 0
                
                if (ovlp_frac1 > max_overlap && ovlp_frac2 > max_overlap)
                    continue
                
                cov1 = (len1 / t_len) * 100
                cov2 = (len2 / t_len) * 100
                
                ref_start = (t_start1 < t_start2) ? t_start1 : t_start2
                ref_end = (t_end1 > t_end2) ? t_end1 : t_end2
                total_cov = ((ref_end - ref_start + 1) / t_len) * 100
                
                if (t_start1 < t_start2) {
                    gap = t_start2 - t_end1 - 1
                    ref_order = "1->2"
                } else {
                    gap = t_start1 - t_end2 - 1
                    ref_order = "2->1"
                }
                
                avg_identity = (pident1 + pident2) / 2
                
                if (total_cov >= 70 && avg_identity >= 50)
                    confidence = "HIGH"
                else if (total_cov >= 50 && avg_identity >= 40)
                    confidence = "MEDIUM"
                else
                    confidence = "LOW"
                
                printf "%s\t%s\t%d-%d\t%d-%d\t%s\t%s\t%d\t%s\t%d\t%.1f\t%.1f\t%.1f\t%.2e\t%.2e\t%d\t%d\t%d\t%d\t%.1f\t%.1f\t%.1f\t%d\t%d\t%s\t%s\n",
                    gene1, gene2, s1, e1, s2, e2, strand1, strand2, distance,
                    t, t_len, pident1, pident2, avg_identity, evalue1, evalue2,
                    t_start1, t_end1, t_start2, t_end2,
                    cov1, cov2, total_cov, gap, overlap, ref_order, confidence
                
                frameshifts++
            }
        }
    }
    
    print "Found", frameshifts, "frameshifts" > "/dev/stderr"
}
' >> "$OUTPUT"

echo "Done!" >&2

