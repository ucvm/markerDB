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
tax_file = snakemake@input$info

rdp_file = snakemake@output$rdp
mothur_file = snakemake@output$mothur
mothur_tax = snakemake@output$mothur_tax
dada2_file = snakemake@output$dada2


# Get sequences and taxonomy info -----------------------------------------

seqs = readDNAStringSet(seq_file)
tax = read_tsv(tax_file)

seq_df = data_frame(
  accn = names(seqs),
  seq = as.character(seqs)
)

db = seq_df %>% left_join(tax)


# Write out proper formats ------------------------------------------------

# --- dada2
dada2 = db %>%
  mutate_all(~str_replace_na(.x)) %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  mutate(tax = str_c(superkingdom, kingdom, phylum, class, order,
                     family, genus, species, sep = ";")) %>%
  pull(tax)

writeXStringSet(setNames(seqs, dada2), dada2_file)

# --- rdp
rdp = db  %>%
  mutate_all(~str_replace_na(.x)) %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  mutate(tax = str_c("Root", superkingdom, kingdom, phylum, class, order,
                     family, genus, species, sep = ";")) %>%
  select(accn, tax) %>%
  mutate(tax = str_c(accn, tax, sep = " ")) %>%
  pull(tax)

writeXStringSet(setNames(seqs, rdp), rdp_file)

# --- mothur
mothur = db  %>%
  mutate_all(~str_replace_na(.x)) %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  mutate(tax = str_c(superkingdom, kingdom, phylum, class, order,
                     family, genus, species, sep = ";")) %>%
  mutate(tax = str_c(tax, ";")) %>%
  select(accn, tax)

writeXStringSet(seqs, mothur_file)
write_tsv(mothur, mothur_tax, col_names = FALSE)





