#!/usr/bin/env Rscript

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to fix ranges when strand is negative
fix_ranges <- function(data) {
    # For each row where strand is "-", swap Start and End if Start > End
    for(i in 1:nrow(data)) {
        if(data$V7[i] == "-" && as.numeric(data$V4[i]) > as.numeric(data$V5[i])) {
            temp <- data$V4[i]
            data$V4[i] <- data$V5[i]
            data$V5[i] <- temp
        }
    }
    return(data)
}

# Function to clean Glimmer file
clean_glimmer <- function(file) {
    data <- read.delim(file, header=FALSE, stringsAsFactors=FALSE)
    # Change format to match required output
    data$V1 <- "contig1_loki"  # Change SeqID
    data$V2 <- "glimmer3"      # Change Source
    data$V6 <- "0"             # Change Score to 0
    data$V8 <- "0"             # Change Phase to 0
    # Fix ranges for negative strand
    data <- fix_ranges(data)
    write.table(data, "cleaned_glimmer.gff", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# Function to clean GeneMark file
clean_genemark <- function(file) {
    lines <- readLines(file)
    # Remove header and footer lines
    data_lines <- lines[!grepl("^#|^##", lines)]
    data_lines <- data_lines[data_lines != ""]
    # Parse the data
    data <- read.delim(text=data_lines, header=FALSE, stringsAsFactors=FALSE)
    data$V1 <- "contig1_loki"  # Change SeqID
    # Fix ranges for negative strand
    data <- fix_ranges(data)
    write.table(data, "cleaned_genemark.gff", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# Function to clean Prodigal file
clean_prodigal <- function(file) {
    lines <- readLines(file)
    output <- character()

    for(line in lines) {
        if(grepl("^>", line)) {
            parts <- strsplit(line, " # ")[[1]]
            id <- gsub("^>", "", parts[1])
            start <- as.numeric(parts[2])
            end <- as.numeric(parts[3])
            strand <- ifelse(parts[4] == "-1", "-", "+")
            
            # Fix ranges if strand is negative
            if(strand == "-" && start > end) {
                temp <- start
                start <- end
                end <- temp
            }
            
            attrs <- gsub(".*ID=", "ID=", parts[5])

            # Prodigal only outputs CDS, so all entries are CDS
            new_line <- paste("contig1_loki",
                            "Prodigal",
                            "CDS",
                            start,
                            end,
                            "0",
                            strand,
                            "0",
                            attrs,
                            sep="\t")
            output <- c(output, new_line)
        }
    }
    writeLines(output, "cleaned_prodigal.gff")
}

# Function to clean Prokka file
clean_prokka <- function(file) {    
    lines <- readLines(file, warn = FALSE)    
    data_lines <- lines[!grepl("^##|^>", lines)]
    data_lines <- data_lines[grepl("\\t", data_lines)]
    
    if(length(data_lines) == 0) {
        stop("ERROR: No valid GFF data lines found in PROKKA file!")
    }
    
    data <- tryCatch({
        read.delim(text=data_lines, header=FALSE, stringsAsFactors=FALSE, comment.char="")
    }, error = function(e) {
        cat("ERROR reading PROKKA data:\n")
        cat("First few lines:\n")
        cat(head(data_lines, 10), sep="\n")
        stop(e)
    })
    
    # Change identifiers
    data$V1 <- "contig1_loki"
    data$V2 <- "Prodigal:2.6"
    
    # Keep ALL feature types (CDS, tRNA, rRNA, etc.)
    # No filtering by feature type
    
    # Fix ranges for negative strand
    data <- fix_ranges(data)
    
    write.table(data, "cleaned_prokka.gff", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# Function to clean RAST file
clean_rast <- function(file) {
    data <- read.delim(file, header=FALSE, skip=1, stringsAsFactors=FALSE)    
    data$V1 <- "contig1_loki"
    # Replace "." with "0" in Score column (V6)
    data$V6[data$V6 == "."] <- "0"
    
    # Keep ALL feature types (CDS, tRNA, rRNA, etc.)
    # No filtering by feature type
    
    # Fix ranges for negative strand
    data <- fix_ranges(data)
    
    # Remove duplicates (RAST has triplicated tRNA/rRNA entries)
    data <- unique(data)
    
    write.table(data, "cleaned_rast.gff", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# Process files according to argument
if(length(args) >= 2) {
    file <- args[2]
    switch(args[1],
           "glimmer" = clean_glimmer(file),
           "genemark" = clean_genemark(file),
           "prodigal" = clean_prodigal(file),
           "prokka" = clean_prokka(file),
           "rast" = clean_rast(file))
} else {
    cat("Usage: Rscript clean_gff.R <tool> <input_file>\n")
    cat("Tools: glimmer, genemark, prodigal, prokka, rast\n")
}

