# ----------------------------------
#  Search for ITS2 sequences
# ----------------------------------

library(rentrez)
library(taxize)
library(glue)
library(purrr)
library(stringr)
library(dplyr)
library(readr)
library(tidyr)


organism = snakemake@config$organism


# search nucleotide using provided organism -------------------------------

term = glue("{organism}[ORGN] (ITS2 OR internal transcribed spacer 2)")
search = entrez_search(db = "nucleotide", term = term, use_history = TRUE)

message("Found ", search$count, " hits")

# Get summaries -----------------------------------------------------------
# needed to we can add taxonomy

message("Getting summaries...")
summaries = seq(0, search$count, by = 200) %>%
  map(~entrez_summary(db = "nucleotide",
                      web_history = search$web_history,
                      retmax = 200,
                      retstart = .x)) %>%
  flatten()

# turn search summaries into a dataframe
sum_df = map_dfr(summaries, ~data_frame(
  uid = as.character(.x[["uid"]]),
  taxid = as.character(.x[["taxid"]])))

# get accension numbers
ids = split(sum_df$uid, ceiling(seq_along(sum_df$uid) / 200))
accns = map(ids, ~entrez_fetch("nucleotide", .x, rettype = "acc"))
accns = str_c(accns, collapse = "\n")
accns = str_split(accns, "\n")[[1]]
accns = discard(accns, ~.x == "")

if (length(accns) != nrow(sum_df)) {
  stop("Number of accession ids retrieved does not match the number of search hits")
}

sum_df = sum_df %>% mutate(accn = accns)


# Get taxonomic classifications -------------------------------------------

message("Retrieving taxonomic classifications...")

taxids = map(summaries, "taxid") %>% unique()
ranks = c("superkingdom", "kingdom", "phylum", "class", "order",
          "family", "genus", "species" )

tax = classification(taxids, db = "ncbi") %>% map(as_data_frame)

tax = tax %>%
  bind_rows(.id = "taxid") %>%
  select(-id) %>%
  filter(rank %in% ranks) %>%
  spread(rank, name) %>%
  select(taxid, ranks)

# add in the taxonomies and write it out
sum_df = sum_df %>% left_join(tax)
message("Writing sequence info")
write_tsv(sum_df, snakemake@output$seqinfo)

# Get fasta sequences -----------------------------------------------------

# ids = split(sum_df$uid, ceiling(seq_along(ids) / 150))
# seqs = map(ids, ~entrez_fetch("nucleotide", id = .x, rettype = "fasta", retmode = "text"))

seqs = seq(0, search$count, by = 50) %>%
  map(~entrez_fetch(db = "nucleotide",
                    web_history = search$web_history,
                    retmax = 50, retstart = .x,
                    rettype = "fasta", retmode = "text")) %>%
  flatten()


message("Writing fasta of hits")
seqs %>%
  walk(~write_file(.x, snakemake@output$seqs, append = TRUE))


message("All done. Goodbye")

