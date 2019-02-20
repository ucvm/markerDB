# ----------------------------------
#  Write out final database
# ----------------------------------

# Writes out database in necessary formats with and without taxonomy info

library(Biostrings)
library(readr)
library(dplyr)
library(stringr)


# Inputs and Outputs ------------------------------------------------------

seq_file = snakemake@input$seqs
tax_file = snakemake@input$taxa

rdp_file = snakemake@output$rdp
mothur_file = snakemake@output$mothur
mothur_tax = snakemake@output$mothur_tax
dada2_file = snakemake@output$dada2


# Get sequences and taxonomy info -----------------------------------------

seqs = Biostrings::readDNAStringSet(seq_file)
tax = readr::read_tsv(tax_file)

seq_df = dplyr::tibble(
  accn = names(seqs),
  seq = as.character(seqs)
)

db = seq_df %>% dplyr::left_join(tax)


# Write out proper formats ------------------------------------------------

# --- dada2
dada2 = db %>%
  dplyr::mutate_all(~stringr::str_replace_na(.x)) %>%
  dplyr::mutate(species = stringr::str_replace_all(species, " ", "_")) %>%
  dplyr::mutate(tax = stringr::str_c(superkingdom, kingdom, phylum, class, order,
  													family, genus, species, sep = ";")) %>%
  dplyr::pull(tax)

Biostrings::writeXStringSet(setNames(seqs, dada2), dada2_file)

# --- rdp
rdp = db  %>%
  dplyr::mutate_all(~stringr::str_replace_na(.x)) %>%
  dplyr::mutate(species = stringr::str_replace_all(species, " ", "_")) %>%
  dplyr::mutate(tax = stringr::str_c("Root", superkingdom, kingdom, phylum, class, order,
  													family, genus, species, sep = ";")) %>%
  dplyr::select(accn, tax) %>%
  dplyr::mutate(tax = stringr::str_c(accn, tax, sep = " ")) %>%
  dplyr::pull(tax)

Biostrings::writeXStringSet(setNames(seqs, rdp), rdp_file)

# --- mothur
mothur = db  %>%
  dplyr::mutate_all(~stringr::str_replace_na(.x)) %>%
  dplyr::mutate(species = stringr::str_replace_all(species, " ", "_")) %>%
  dplyr::mutate(tax = stringr::str_c(superkingdom, kingdom, phylum, class, order,
  													family, genus, species, sep = ";")) %>%
  dplyr::mutate(tax = stringr::str_c(tax, ";")) %>%
  dplyr::select(accn, tax)

Biostrings::writeXStringSet(seqs, mothur_file)
readr::write_tsv(mothur, mothur_tax, col_names = FALSE)





