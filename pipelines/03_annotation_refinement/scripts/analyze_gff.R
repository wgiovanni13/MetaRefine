#!/usr/bin/env Rscript

.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))

library(ggplot2)
library(plotly)
library(dplyr)
library(plotly)
library(data.table)


gms2 <- read.delim("cleaned_genemark.gff", header = FALSE, col.names = c("SeqID", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"))
glimmer <- read.delim("cleaned_glimmer.gff", header = FALSE, col.names = c("SeqID", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"))
prodigal <- read.delim("cleaned_prodigal.gff", header = FALSE, col.names = c("SeqID", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"))
prokka <- read.delim("cleaned_prokka.gff", header = FALSE, col.names = c("SeqID", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"))
rast <- read.delim("cleaned_rast.gff", header = FALSE, col.names = c("SeqID", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"))


############## Calculate overlap percentage #################

# Function to verify if two ranges has overlap

check_overlap <- function(start1, end1, start2, end2) {
return(!(end1 < start2 || end2 < start1)) # True si hay overlap
}

# Function to compare dataframes

compare_ranges <- function(df_list) {
results <- list() # Almacenar resultados

# Create a list of dataframes combinations

combinations <- expand.grid(1:length(df_list), 1:length(df_list))
combinations <- combinations[combinations$Var1 != combinations$Var2, ]

# Function to process each combination

process_combination <- function(i, j) {
df1 <- df_list[[i]]
df2 <- df_list[[j]]

# Compare each row from df1 against all rows in df2

result <- rbindlist(lapply(1:nrow(df1), function(k) {
overlaps <- sapply(1:nrow(df2), function(l) {
if (check_overlap(df1$Start[k], df1$End[k], df2$Start[l], df2$End[l])) {
return(data.table(
DF1_Source = i,
DF2_Source = j,
SeqID_DF1 = df1$SeqID[k],
Start_DF1 = df1$Start[k],
End_DF1 = df1$End[k],
Strand_DF1 = df1$Strand[k],
Attributes_DF1 = df1$Attributes[k],
SeqID_DF2 = df2$SeqID[l],
Start_DF2 = df2$Start[l],
End_DF2 = df2$End[l],
Strand_DF2 = df2$Strand[l],
Attributes_DF2 = df2$Attributes[l]
))
} else {
return(NULL)
}
})
return(rbindlist(overlaps, use.names = TRUE, fill = TRUE))
}), use.names = TRUE, fill = TRUE)

return(result)
}

# Process all combinations

for (row in 1:nrow(combinations)) {
i <- combinations[row, 1]
j <- combinations[row, 2]
results[[length(results) + 1]] <- process_combination(i, j)
}

# Convert all list of results in a single dataframe

return(rbindlist(results, use.names = TRUE, fill = TRUE))
}

# List of dataframes

df_list <- list(gms2, glimmer, prodigal, prokka, rast)

# Execute comparison

overlap_results <- compare_ranges(df_list)

dfs <- list(
gms2 = as.data.table(gms2),
glimmer = as.data.table(glimmer),
prodigal = as.data.table(prodigal),
prokka = as.data.table(prokka),
rast = as.data.table(rast)
)

# Function to look for bidirectional overlaps with percentage of overlapping

find_bidirectional_overlaps <- function(dt1, dt2, source1, source2) {

# Establish keys for efficient unions

setkey(dt1, Start, End)
setkey(dt2, Start, End)

# Overlap: dt1 against dt2

overlaps_1_2 <- foverlaps(
dt1[, .(Start, End, SeqID, Strand, Attributes)],
dt2[, .(Start, End, SeqID, Strand, Attributes)],
by.x = c("Start", "End"),
type = "any", # Look for any overlap
nomatch = 0L # Exclude rows without coincidences
)
overlaps_1_2[, `:=`(
DF1_Source = source1,
DF2_Source = source2,
Overlap_Length = pmin(i.End, End) - pmax(i.Start, Start) + 1,
DF1_Length = i.End - i.Start + 1,
DF2_Length = End - Start + 1
)]
overlaps_1_2[, `:=`(
Overlap_Percentage_DF1 = (Overlap_Length / DF1_Length) * 100,
Overlap_Percentage_DF2 = (Overlap_Length / DF2_Length) * 100
)]

# Overlap: dt2 against dt1

overlaps_2_1 <- foverlaps(
dt2[, .(Start, End, SeqID, Strand, Attributes)],
dt1[, .(Start, End, SeqID, Strand, Attributes)],
by.x = c("Start", "End"),
type = "any",
nomatch = 0L
)
overlaps_2_1[, `:=`(
DF1_Source = source2,
DF2_Source = source1,
Overlap_Length = pmin(i.End, End) - pmax(i.Start, Start) + 1,
DF1_Length = i.End - i.Start + 1,
DF2_Length = End - Start + 1
)]
overlaps_2_1[, `:=`(
Overlap_Percentage_DF1 = (Overlap_Length / DF1_Length) * 100,
Overlap_Percentage_DF2 = (Overlap_Length / DF2_Length) * 100
)]

# Combine results

combined <- rbind(
overlaps_1_2[, .(
DF1_Source, DF2_Source,
SeqID_DF1 = SeqID, Start_DF1 = i.Start, End_DF1 = i.End, Strand_DF1 = i.Strand, Attributes_DF1 = i.Attributes,
SeqID_DF2 = SeqID, Start_DF2 = Start, End_DF2 = End, Strand_DF2 = Strand, Attributes_DF2 = Attributes,
Overlap_Length, Overlap_Percentage_DF1, Overlap_Percentage_DF2
)],
overlaps_2_1[, .(
DF1_Source, DF2_Source,
SeqID_DF1 = SeqID, Start_DF1 = Start, End_DF1 = End, Strand_DF1 = Strand, Attributes_DF1 = Attributes,
SeqID_DF2 = SeqID, Start_DF2 = i.Start, End_DF2 = i.End, Strand_DF2 = i.Strand, Attributes_DF2 = i.Attributes,
Overlap_Length, Overlap_Percentage_DF1, Overlap_Percentage_DF2
)]
)

return(combined)
}

# Do comparisons between all pairs from data.table

results <- rbindlist(lapply(names(dfs), function(df1_name) {
rbindlist(lapply(names(dfs), function(df2_name) {
if (df1_name != df2_name) {
find_bidirectional_overlaps(dfs[[df1_name]], dfs[[df2_name]], df1_name, df2_name)
} else {
NULL
}
}))
}))


results_unique <- results[
, `:=`(
Unique_ID = paste(
pmin(DF1_Source, DF2_Source),
pmax(DF1_Source, DF2_Source),
pmin(SeqID_DF1, SeqID_DF2),
pmax(SeqID_DF1, SeqID_DF2),
pmin(Start_DF1, Start_DF2),
pmax(Start_DF1, Start_DF2),
sep = "_"
)
)
][!duplicated(Unique_ID)]

# Sequences that don't have any overlap with others

missing_sequences <- list()

# Iterate over the names of the dataframes

for (tool_name in names(dfs)) {

# Select columns Attributes_DF1 Attributes_DF2 where the tools appear

relevant_attributes <- c(
results_unique[DF1_Source == tool_name, Attributes_DF1],
results_unique[DF2_Source == tool_name, Attributes_DF2]
)

# Obtain unique sequences

relevant_attributes <- unique(relevant_attributes)

# Identify sequences that are in the original file but not in the results_unique

original_attributes <- unique(dfs[[tool_name]]$Attributes)
missing_attributes <- setdiff(original_attributes, relevant_attributes)

missing_sequences[[tool_name]] <- dfs[[tool_name]][Attributes %in% missing_attributes]
}

for (tool_name in names(missing_sequences)) {
fwrite(
missing_sequences[[tool_name]],
paste0("missing_sequences_", tool_name, ".csv")
)
}

gms2_missing_sequences <- missing_sequences[["gms2"]]
glimmer_missing_sequences <- missing_sequences[["glimmer"]]
prodigal_missing_sequences <- missing_sequences[["prodigal"]]
prokka_missing_sequences <- missing_sequences[["prokka"]]
rast_missing_sequences <- missing_sequences[["rast"]]

# Add a column with 0% of overlap to this new dataframes with unique sequences

gms2_missing_sequences[, Overlap := 0.00]
glimmer_missing_sequences[, Overlap := 0.00]
prodigal_missing_sequences[, Overlap := 0.00]
prokka_missing_sequences[, Overlap := 0.00]
rast_missing_sequences[, Overlap := 0.00]


########## Categorizing overlapping sequences in groups ##########

# Step 1: Categorizing all the sequences that are overlapping

ranges_combined <- rbindlist(list(
results_unique[, .(Start = Start_DF1, End = End_DF1, Tool = DF1_Source)],
results_unique[, .(Start = Start_DF2, End = End_DF2, Tool = DF2_Source)]
))

# Step 2: Groupping in unique ranges and consolidate the tools belonging on them

range_summary <- ranges_combined[, .(
Tools = sort(unique(Tool)), # Ordered list of unique tools
Tool_Count = uniqueN(Tool) # Nº of unique tools
), by = .(Start, End)]

# Step 3: Creating a list of tools as strings

range_summary[, Tools_List := paste(Tools, collapse = ", ")]

# Step 4: Filter ranges according the quantity of tool numbers

ranges_in_5_tools <- range_summary[Tool_Count == 5, .(Start, End, Tools_List)]
ranges_in_4_tools <- range_summary[Tool_Count == 4, .(Start, End, Tools_List)]
ranges_in_3_tools <- range_summary[Tool_Count == 3, .(Start, End, Tools_List)]
ranges_in_2_tools <- range_summary[Tool_Count == 2, .(Start, End, Tools_List)]
ranges_in_1_tool <- range_summary[Tool_Count == 1, .(Start, End, Tools_List)]


clean_dataframe <- function(df) {

# Delete duplicates based on Start and End positions

df <- unique(df, by = c("Start", "End"))

# Assure Tools_List has unique and ordered tools

df[, Tools_List := paste(sort(unique(unlist(strsplit(Tools_List, ", ")))), collapse = ", ")]
return(df)
}

# Apply the "cleaning to each filtered tools dataframe

ranges_in_5_tools <- clean_dataframe(ranges_in_5_tools)
ranges_in_4_tools <- clean_dataframe(ranges_in_4_tools)
ranges_in_3_tools <- clean_dataframe(ranges_in_3_tools)
ranges_in_2_tools <- clean_dataframe(ranges_in_2_tools)
ranges_in_1_tool <- clean_dataframe(ranges_in_1_tool)

ranges_in_5_tools
ranges_in_4_tools
ranges_in_3_tools
ranges_in_2_tools
ranges_in_1_tool

# Function to identify the tool based on the range

identify_tool <- function(start, end, dfs) {
tool_names <- names(dfs)

for (tool in tool_names) {
tool_df <- dfs[[tool]]

# Looking for exact coincidences of the ranges in the dataframe of the tools

if (any(tool_df$Start == start & tool_df$End == end)) {
return(tool) # Return the name of the tool
}
}
return(NA) # If there are no tools assign to the range, it returns NA
}

# Update Tools_List column in ranges_in_1_tool

ranges_in_1_tool[, Tools_List := sapply(1:nrow(ranges_in_1_tool), function(idx) {
identify_tool(ranges_in_1_tool$Start[idx], ranges_in_1_tool$End[idx], dfs)
})]

# Function to identify the tools belonging in "ranges_in_2_tools" "ranges_in_3_tools" "ranges_in_4_tools"

update_tools_from_results <- function(df, results_df) {
df[, Tools_List := sapply(1:nrow(df), function(idx) {
start <- df$Start[idx]
end <- df$End[idx]

# Look for tools in results_unique with coincidents ranges

tools <- unique(c(
results_df[(Start_DF1 == start & End_DF1 == end), DF1_Source],
results_df[(Start_DF2 == start & End_DF2 == end), DF2_Source]
))

# Return tools concatenated

paste(tools, collapse = ", ")
})]
return(df)
}

# Apply function to the dataframes

ranges_in_5_tools <- update_tools_from_results(ranges_in_5_tools, results_unique)
ranges_in_4_tools <- update_tools_from_results(ranges_in_4_tools, results_unique)
ranges_in_3_tools <- update_tools_from_results(ranges_in_3_tools, results_unique)
ranges_in_2_tools <- update_tools_from_results(ranges_in_2_tools, results_unique)

# Order all the 5 dataframes from lower to greater value based on first column

range_dfs <- list(
ranges_in_5_tools,
ranges_in_4_tools,
ranges_in_3_tools,
ranges_in_2_tools,
ranges_in_1_tool
)

# Order each dataframe by first column (start position)

range_dfs <- lapply(range_dfs, function(df) {
setorder(df, Start)
return(df)
})

# Rename dataframes to keep original variables and to not have any mess

ranges_in_5_tools <- range_dfs[[1]]
ranges_in_4_tools <- range_dfs[[2]]
ranges_in_3_tools <- range_dfs[[3]]
ranges_in_2_tools <- range_dfs[[4]]
ranges_in_1_tool <- range_dfs[[5]]

# Adding strand and feature for ranges_in_1_tool, ranges_in_2_tools, ranges_in_3_tools, ranges_in_4_tools, ranges_in_5_tools

tools_dataframes <- list(
  rast = rast,
  prodigal = prodigal,
  prokka = prokka,
  glimmer = glimmer,
  gms2 = gms2
)

# Make sure columns have correct names for tools 

for (tool_name in names(tools_dataframes)) {
  tool_df <- tools_dataframes[[tool_name]]
  
  if (!all(c("Start", "End", "Strand", "Type") %in% colnames(tool_df))) {
    stop(paste("Dataframe of the tool", tool_name, "does not have necessary columns: Start, End, Strand, Type"))
  }
}

# Function to process dataframes with ranges

process_ranges <- function(ranges_df) {

  # Create columns needed for all the ranges Strand_List and Feature_List

  ranges_df$Strand_List <- NA
  ranges_df$Feature_List <- NA
  
  # A for look for each row cada fila del dataframe de rangos

  for (i in seq_len(nrow(ranges_df))) {
    start <- as.numeric(ranges_df[i, 1])  # Column 1: Start of the range
    end <- as.numeric(ranges_df[i, 2])    # Column 2: End of the range
    tools <- as.character(ranges_df[i, 3])  # Column 3: Tools
    
    # Separate tools by coma

    tool_list <- unlist(strsplit(tools, split = ","))
    
    # Take the first tool

    primary_tool <- trimws(tool_list[1])
    
    # Check in the correspond dataframe

    if (primary_tool %in% names(tools_dataframes)) {
      tool_df <- tools_dataframes[[primary_tool]]
      
      # Filter rows in the dataframe of the tool 

      match <- tool_df %>%
        filter(Start <= start & End >= end)  # Column in the ranges of the dataframe of the tool
      
      if (nrow(match) > 0) {

        # Take value in the columns strand and type to add in the new two columns Strand_List and Feature_List

        ranges_df$Strand_List[i] <- match$Strand[1]
        ranges_df$Feature_List[i] <- match$Type[1]
      }
    }
  }
  
  return(ranges_df)
}

# Apply function to each dataframe of ranges

ranges_in_1_tool <- process_ranges(ranges_in_1_tool)
ranges_in_2_tools <- process_ranges(ranges_in_2_tools)
ranges_in_3_tools <- process_ranges(ranges_in_3_tools)
ranges_in_4_tools <- process_ranges(ranges_in_4_tools)
ranges_in_5_tools <- process_ranges(ranges_in_5_tools)

############## Calculate overlap between groups of ranges_in_1_tool, ranges_in_2_tools, ranges_in_3_tools, ranges_in_4_tools, ranges_in_5_tools #################

compare_ranges_with_asymmetric_overlaps <- function(df1, df2, include_self_comparison = FALSE, include_all_from_df1 = FALSE) {
results <- data.table()

for (i in 1:nrow(df1)) {
range1_start <- df1$Start[i]
range1_end <- df1$End[i]
strand1 <- df1$Strand_List[i]
feature1 <- df1$Feature_List[i]
tools1 <- df1$Tools_List[i]

overlaps <- data.table()

for (j in 1:nrow(df2)) {
range2_start <- df2$Start[j]
range2_end <- df2$End[j]
strand2 <- df2$Strand_List[j]
feature2 <- df2$Feature_List[j]
tools2 <- df2$Tools_List[j]

# Calculate overlapping

overlap_start <- max(range1_start, range2_start)
overlap_end <- min(range1_end, range2_end)

if (overlap_start <= overlap_end) {
overlap_length <- overlap_end - overlap_start + 1
range1_length <- range1_end - range1_start + 1
range2_length <- range2_end - range2_start + 1

# Percentage of overlaps from both sides

overlap_percentage_df1 <- 100 * overlap_length / range1_length
overlap_percentage_df2 <- 100 * overlap_length / range2_length

# Filter combinations

if (!include_self_comparison || i != j || overlap_percentage_df1 > 0) {
overlaps <- rbind(overlaps, data.table(
Start_DF1 = range1_start,
End_DF1 = range1_end,
Start_DF2 = range2_start,
End_DF2 = range2_end,
Overlap_Percentage_DF1 = overlap_percentage_df1,
Overlap_Percentage_DF2 = overlap_percentage_df2,
Strand_DF1 = strand1,
Strand_DF2 = strand2,
Feature_DF1 = feature1,
Feature_DF2 = feature2,
Tools_DF1 = tools1,
Tools_DF2 = tools2
))
}
}
}

# If there are no overlaps and this requires including all sequences from df1

if (include_all_from_df1 && nrow(overlaps) == 0) {
overlaps <- data.table(
Start_DF1 = range1_start,
End_DF1 = range1_end,
Start_DF2 = NA,
End_DF2 = NA,
Overlap_Percentage_DF1 = 0,
Overlap_Percentage_DF2 = NA,
Strand_DF1 = strand1,
Strand_DF2 = NA,
Feature_DF1 = feature1,
Feature_DF2 = NA,
Tools_DF1 = tools1,
Tools_DF2 = NA
)
}

# Add results

results <- rbind(results, overlaps)
}

return(results)
}

# Comparisons among all dataframes

all_dataframes <- list(ranges_in_5_tools, ranges_in_4_tools, ranges_in_3_tools, ranges_in_2_tools, ranges_in_1_tool)

comparison_results <- list()

for (i in 1:length(all_dataframes)) {
for (j in 1:length(all_dataframes)) {
include_self_comparison <- (i == j && i == length(all_dataframes))
include_all_from_df1 <- (i == 1 && j == 1)
comparison <- compare_ranges_with_asymmetric_overlaps(all_dataframes[[i]], all_dataframes[[j]],
include_self_comparison, include_all_from_df1)
if (nrow(comparison) > 0) {
comparison_results <- append(comparison_results, list(comparison))
}
}
}

# Combine all results

final_overlap_results <- rbindlist(comparison_results, fill = TRUE)

final_overlap_results <- final_overlap_results %>%
arrange(Start_DF1, -End_DF1)


# Filter unique ranges in the first two columns


# Change information for all the sequences founded by 5 tools

full_tools <- "gms2, glimmer, prodigal, prokka, rast"

# Modify files from dataframe final_overlap_results

final_overlap_results[

# Condition: equal ranges and both columns from tools which coincided with all the tools

Start_DF1 == Start_DF2 & End_DF1 == End_DF2 &
Tools_DF1 == full_tools & Tools_DF2 == full_tools,

# Modify all rows selected

`:=`(
Start_DF2 = NA,
End_DF2 = NA,
Overlap_Percentage_DF1 = 100,
Overlap_Percentage_DF2 = 0,
Strand_DF2 = NA,
Feature_DF2 = NA,
Tools_DF2 = NA
)
]


# ============================================================================
# REMOVE DUPLICATES (SMART VERSION)
# ============================================================================

cat("\n================================================================================\n")
cat("REMOVING DUPLICATES (SMART VERSION)\n")
cat("================================================================================\n")

# Set as a Data Table
setDT(final_overlap_results)

# Clean up any existing auxiliary columns from previous runs
if ("n_combinations" %in% names(final_overlap_results)) {
  final_overlap_results[, n_combinations := NULL]
}
if ("DF2_combo" %in% names(final_overlap_results)) {
  final_overlap_results[, DF2_combo := NULL]
}
if ("is_potential_duplicate" %in% names(final_overlap_results)) {
  final_overlap_results[, is_potential_duplicate := NULL]
}
if ("is_real_duplicate" %in% names(final_overlap_results)) {
  final_overlap_results[, is_real_duplicate := NULL]
}

# Step 1: Identify rows where all columns are identical (potential duplicates)
final_overlap_results[, is_potential_duplicate :=
                        (Start_DF1 == Start_DF2 & End_DF1 == End_DF2) &
                        (Overlap_Percentage_DF1 == Overlap_Percentage_DF2) &
                        (Strand_DF1 == Strand_DF2) &
                        (Feature_DF1 == Feature_DF2) &
                        (Tools_DF1 == Tools_DF2)
]

cat("Potential duplicates identified:", sum(final_overlap_results$is_potential_duplicate, na.rm = TRUE), "\n")

# Step 2: For each range (Start_DF1, End_DF1, Tools_DF1), count how many 
# DIFFERENT DF2 combinations exist

# Create a unique identifier for each DF2 combination (handling NAs)
final_overlap_results[, DF2_combo := paste(
  ifelse(is.na(Start_DF2), "NA", as.character(Start_DF2)),
  ifelse(is.na(End_DF2), "NA", as.character(End_DF2)),
  ifelse(is.na(Tools_DF2), "NA", as.character(Tools_DF2)),
  sep = "_"
)]

# For each range, count unique DF2 combinations
range_info <- final_overlap_results[, .(
  n_combinations = uniqueN(DF2_combo)
), by = .(Start_DF1, End_DF1, Tools_DF1)]

# Merge back to original data
final_overlap_results <- merge(
  final_overlap_results,
  range_info,
  by = c("Start_DF1", "End_DF1", "Tools_DF1"),
  all.x = TRUE
)

# Mark as real duplicate if:
# 1. It's a potential duplicate (all columns match)
# 2. AND there are multiple different DF2 combinations for this range
final_overlap_results[, is_real_duplicate := 
                        is_potential_duplicate & n_combinations > 1
]

real_duplicates <- sum(final_overlap_results$is_real_duplicate, na.rm = TRUE)
false_positives <- sum(final_overlap_results$is_potential_duplicate & !final_overlap_results$is_real_duplicate, na.rm = TRUE)

cat("Real duplicates to remove:", real_duplicates, "\n")
cat("False positives (will be kept):", false_positives, "\n")

# Step 3: Remove only real duplicates
final_overlap_results_clean <- final_overlap_results[is_real_duplicate == FALSE]

# Clean up auxiliary columns
final_overlap_results_clean[, `:=`(
  is_potential_duplicate = NULL, 
  is_real_duplicate = NULL,
  DF2_combo = NULL,
  n_combinations = NULL
)]

cat("Rows after removing duplicates:", nrow(final_overlap_results_clean), "\n")


# ============================================================================
# REMOVE INVERSE ROWS
# ============================================================================

cat("\n================================================================================\n")
cat("REMOVING INVERSE ROWS\n")
cat("================================================================================\n")

setDT(final_overlap_results_clean)

# Clean up any existing auxiliary columns from previous runs
if ("inversion_key" %in% names(final_overlap_results_clean)) {
  final_overlap_results_clean[, inversion_key := NULL]
}
if ("is_inverted_pair" %in% names(final_overlap_results_clean)) {
  final_overlap_results_clean[, is_inverted_pair := NULL]
}
if ("to_remove" %in% names(final_overlap_results_clean)) {
  final_overlap_results_clean[, to_remove := NULL]
}

before_inversions <- nrow(final_overlap_results_clean)

# Create a unique key that identifies potentially inverted rows
# This key will be the same for both A->B and B->A
final_overlap_results_clean[, inversion_key := paste(
  pmin(Start_DF1, Start_DF2, na.rm = TRUE),   # Min of Starts
  pmax(Start_DF1, Start_DF2, na.rm = TRUE),   # Max of Starts
  pmin(End_DF1, End_DF2, na.rm = TRUE),       # Min of Ends
  pmax(End_DF1, End_DF2, na.rm = TRUE),       # Max of Ends
  pmin(Tools_DF1, Tools_DF2, na.rm = TRUE),   # Min of Tools (alphabetically)
  pmax(Tools_DF1, Tools_DF2, na.rm = TRUE),   # Max of Tools
  sep = "_"
)]

# Identify duplicated rows based on inversion key
# Exclude rows with key "NA_NA_NA_NA_NA_NA" (rows without DF2)
final_overlap_results_clean[, is_inverted_pair := 
                              duplicated(inversion_key) | duplicated(inversion_key, fromLast = TRUE)
]

# Mark to remove only the SECOND occurrence of inverted pairs
final_overlap_results_clean[, to_remove := 
                              is_inverted_pair & 
                              duplicated(inversion_key) & 
                              inversion_key != "NA_NA_NA_NA_NA_NA"
]

inversions_to_remove <- sum(final_overlap_results_clean$to_remove, na.rm = TRUE)
cat("Inverse rows to remove:", inversions_to_remove, "\n")

# Filter rows, removing those with `to_remove == TRUE`
final_overlap_results_clean <- final_overlap_results_clean[to_remove == FALSE]

# Remove auxiliary columns
final_overlap_results_clean[, `:=`(
  inversion_key = NULL, 
  is_inverted_pair = NULL, 
  to_remove = NULL
)]

cat("Rows after removing inversions:", nrow(final_overlap_results_clean), "\n")
cat("Total rows removed:", before_inversions - nrow(final_overlap_results_clean), "\n")

# ============================================================================
# ADD MISSING SEQUENCES TO FINAL DATAFRAME
# ============================================================================

cat("\n================================================================================\n")
cat("ADDING MISSING SEQUENCES TO FINAL DATAFRAME\n")
cat("================================================================================\n")

# Combine all missing sequences into one data.table
missing_sequences_combined <- rbindlist(missing_sequences)

cat("Total missing sequences before cleaning:", nrow(missing_sequences_combined), "\n")

# Select only necessary columns
missing_sequences_combined <- missing_sequences_combined[, .(SeqID, Start, End, Source, Strand, Type)]

# Sort by Start and End (descending)
missing_sequences_combined <- missing_sequences_combined %>%
  arrange(Start, -End)

cat("Missing sequences after selecting columns:", nrow(missing_sequences_combined), "\n")

# Standardize tool names in missing_sequences dataframe
setDT(missing_sequences_combined)

# Option 1: Use fifelse() chain (recommended)
missing_sequences_combined[, Source := fifelse(
  Source == "GeneMark.hmm2", "gms2",
  fifelse(Source == "FIG", "rast",
          fifelse(Source == "glimmer3", "glimmer",
                  fifelse(Source %in% c("Prodigal:2.6", "Aragorn:1.2"), "prokka",
                          fifelse(Source == "Prodigal", "prodigal",
                                  Source))))
)]

cat("Tool name standardization complete\n")

# Create expanded version with structure matching final_overlap_results_clean
missing_sequences_expanded <- data.table(
  Start_DF1 = missing_sequences_combined$Start,
  End_DF1 = missing_sequences_combined$End,
  Start_DF2 = NA_integer_,
  End_DF2 = NA_integer_,
  Overlap_Percentage_DF1 = 0.0,
  Overlap_Percentage_DF2 = 0.0,
  Strand_DF1 = missing_sequences_combined$Strand,
  Strand_DF2 = NA_character_,
  Feature_DF1 = missing_sequences_combined$Type,
  Feature_DF2 = NA_character_,
  Tools_DF1 = missing_sequences_combined$Source,
  Tools_DF2 = NA_character_
)

cat("Missing sequences expanded:", nrow(missing_sequences_expanded), "\n")

# Ensure column types match before combining
# This is CRITICAL to avoid type mismatches during rbind

# Check and align column types
cat("\n--- Checking column type compatibility ---\n")

# Get column types from both data.tables
types_clean <- sapply(final_overlap_results_clean, class)
types_missing <- sapply(missing_sequences_expanded, class)

# Print column types for verification
cat("\nColumn types in final_overlap_results_clean:\n")
print(types_clean)

cat("\nColumn types in missing_sequences_expanded:\n")
print(types_missing)

# Align numeric columns (ensure both are numeric)
numeric_cols <- c("Start_DF1", "End_DF1", "Start_DF2", "End_DF2", 
                  "Overlap_Percentage_DF1", "Overlap_Percentage_DF2")

for (col in numeric_cols) {
  if (col %in% names(final_overlap_results_clean)) {
    final_overlap_results_clean[[col]] <- as.numeric(final_overlap_results_clean[[col]])
  }
  if (col %in% names(missing_sequences_expanded)) {
    missing_sequences_expanded[[col]] <- as.numeric(missing_sequences_expanded[[col]])
  }
}

# Align character columns (ensure both are character)
char_cols <- c("Strand_DF1", "Strand_DF2", "Feature_DF1", "Feature_DF2", 
               "Tools_DF1", "Tools_DF2")

for (col in char_cols) {
  if (col %in% names(final_overlap_results_clean)) {
    final_overlap_results_clean[[col]] <- as.character(final_overlap_results_clean[[col]])
  }
  if (col %in% names(missing_sequences_expanded)) {
    missing_sequences_expanded[[col]] <- as.character(missing_sequences_expanded[[col]])
  }
}

cat("\n✓ Column types aligned\n")

# Concatenate final_overlap_results_clean with missing_sequences_expanded
final_overlap_results_combined <- rbind(
  final_overlap_results_clean, 
  missing_sequences_expanded,
  fill = TRUE
)

cat("\n✓ Combined dataframes\n")
cat("  Rows from final_overlap_results_clean:", nrow(final_overlap_results_clean), "\n")
cat("  Rows from missing_sequences_expanded:", nrow(missing_sequences_expanded), "\n")
cat("  Total rows in combined dataframe:", nrow(final_overlap_results_combined), "\n")

# Sort by Start_DF1 and End_DF1 (descending)
setorder(final_overlap_results_combined, Start_DF1, -End_DF1)

cat("\n✓ Dataframe sorted\n")

# ============================================================================
# EXPORT FINAL COMBINED FILE
# ============================================================================

cat("\n================================================================================\n")
cat("EXPORTING FINAL COMBINED FILE\n")
cat("================================================================================\n")

# Ensure column order is correct
column_order <- c(
  "Start_DF1", "End_DF1", "Start_DF2", "End_DF2",
  "Overlap_Percentage_DF1", "Overlap_Percentage_DF2",
  "Strand_DF1", "Strand_DF2",
  "Feature_DF1", "Feature_DF2",
  "Tools_DF1", "Tools_DF2"
)

# Reorder columns if necessary
if (all(column_order %in% names(final_overlap_results_combined))) {
  final_overlap_results_combined <- final_overlap_results_combined[, ..column_order]
  cat("✓ Columns reordered to match expected format\n")
}

# Export to TSV
write.table(
  final_overlap_results_combined, 
  "final_overlap_results_combined.tsv", 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE,
  na = "NA"
)

cat("\n✓ File exported: final_overlap_results_combined.tsv\n")
cat("  Total rows:", nrow(final_overlap_results_combined), "\n")
cat("  Total columns:", ncol(final_overlap_results_combined), "\n")

