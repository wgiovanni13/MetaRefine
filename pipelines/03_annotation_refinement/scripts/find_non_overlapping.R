#!/usr/bin/env Rscript

# Load libraries

library(data.table)

# Import files

prokka_gff <- fread("cleaned_prokka.gff", 
                    sep="\t",
                    col.names=c("SeqID","Source","Type","Start","End","Score","Strand","Phase","Attributes"))

combined_results <- fread("final_overlap_results_combined.tsv", 
                         sep="\t")

# Extracting ranges from prokka_gff

prokka_ranges <- prokka_gff[, .(Start, End)]

# Function to find no overlaps with prokka

no_overlap_any <- function(start_query, end_query) {
  if (is.na(start_query) || is.na(end_query)) return(FALSE)
  
  all(
    start_query > prokka_ranges$End | 
      end_query < prokka_ranges$Start
  )
}

# Apply for both ranges in combined results

combined_results[, No_Overlap_DF1 := mapply(
  no_overlap_any,
  Start_DF1,
  End_DF1
)]

combined_results[, No_Overlap_DF2 := mapply(
  no_overlap_any,
  Start_DF2,
  End_DF2
)]

# Extract non overlapping ranges excluding NAs

non_overlapping_ranges <- rbind(
  combined_results[No_Overlap_DF1 == TRUE, .(Start = Start_DF1, End = End_DF1, Source = "DF1")],
  combined_results[No_Overlap_DF2 == TRUE, .(Start = Start_DF2, End = End_DF2, Source = "DF2")]
)

# Order by start position

setorder(non_overlapping_ranges, Start)

# Delete possible duplicates

non_overlapping_ranges <- unique(non_overlapping_ranges)

# Add length column

non_overlapping_ranges[, Length := End - Start] 
setorder(non_overlapping_ranges, Start, -Length)

# Identify ranges inside others

non_overlapping_ranges[, Contained := FALSE]

for (i in 1:nrow(non_overlapping_ranges)) {
  if (non_overlapping_ranges$Contained[i]) next
  
  current_start <- non_overlapping_ranges$Start[i]
  current_end <- non_overlapping_ranges$End[i]
  
  contained_indices <- which(
    non_overlapping_ranges$Start >= current_start &
      non_overlapping_ranges$End <= current_end &
      !non_overlapping_ranges$Contained &
      seq_len(nrow(non_overlapping_ranges)) > i
  )
  
  if (length(contained_indices) > 0) {
    non_overlapping_ranges[contained_indices, Contained := TRUE]
  }
}

# Filter ranges not contained in others

final_ranges <- non_overlapping_ranges[Contained == FALSE, .(Start, End, Source)]

# Order by start position 

setorder(final_ranges, Start)

# Save results

fwrite(final_ranges, "non_overlapping_ranges.tsv", sep="\t")

