# ----------------------------------
#  Write out final database
# ----------------------------------

# Writes out database in necessary formats with and without taxonomy info

library(Biostrings)
library(readr)
library(dplyr)
library(stringr)
library(purrr)

source("scripts/write_functions.R")

# Inputs and Outputs ------------------------------------------------------

seq_file = snakemake@input$seqs
tax_file = snakemake@input$taxa
outdir = snakemake@output[[1]]

if (exists(snakemake@input$aln)) {
	aln = snakemake@input$aln
}


dir.create(outdir, showWarnings = FALSE)


# Get sequences and taxonomy info -----------------------------------------

seqs = Biostrings::readDNAStringSet(seq_file)
tax = readr::read_tsv(tax_file)

seq_df = dplyr::tibble(
  accn = names(seqs),
  seq = as.character(seqs)
)

db = seq_df %>% dplyr::left_join(tax)


# Write out proper formats ------------------------------------------------

if (exists(aln)) {
	invoke_map(write_functions, list(list(db = db, seqs = seqs, outdir = outdir, align = aln)))
} else {
	invoke_map(write_functions, list(list(db = db, seqs = seqs, outdir = outdir)))
}





