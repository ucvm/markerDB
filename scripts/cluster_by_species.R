# ----------------------------------
#  Write out final database
# ----------------------------------

# Writes out database in necessary formats with and without taxonomy info

library(Biostrings)
library(readr)
library(dplyr)
library(stringr)
library(glue)
library(purrr)


# Inputs and Outputs ------------------------------------------------------

seq_file = snakemake@input$seqs
tax_file = snakemake@input$info

rdp_file = snakemake@output$rdp
mothur_file = snakemake@output$mothur
mothur_tax = snakemake@output$mothur_tax
dada2_file = snakemake@output$dada2

clusterid = snakemake@config$clusterid

# Get sequences and taxonomy info -----------------------------------------

seqs = readDNAStringSet(seq_file)
tax = read_tsv(tax_file)

seq_df = data_frame(
  accn = names(seqs),
  seq = as.character(seqs)
)

db = seq_df %>% left_join(tax)

by_species = db %>% split(db$taxid)

message("Clustering each species at ", round(clusterid * 100, digits = 0), " % identity")

cluster_res = by_species %>% imap(function(t, tax) {
  seq_file = tempfile()
  writeXStringSet(DNAStringSet(setNames(t$seq, t$accn), use.names = TRUE), seq_file)
  cmd = "vsearch"
  outfile = tempfile()
  args = glue("--cluster_fast {seq_file} -id {clusterid} --clusterout_sort --consout - | seqkit head -n1 ")
  system2(cmd, args = args, stdout = outfile, stderr = NULL) 
  seqs = readDNAStringSet(outfile)
  names(seqs) = tax
  return(seqs)
})

clustered_seqs = map(cluster_res, ~.x[[1]])

clustered_seqs = DNAStringSet(clustered_seqs, use.names = TRUE)

tax_only = tax %>% 
  distinct(taxid, superkingdom, kingdom, phylum, class, order, family, genus, species) %>% 
  mutate(taxid = as.character(taxid))

db_clus = data_frame(taxid = names(clustered_seqs)) %>% 
  left_join(tax_only)



# Write out proper formats ------------------------------------------------

message("Writing out clustered data")

# --- dada2
dada2 = db_clus %>%
  mutate_all(~str_replace_na(.x)) %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  mutate(tax = str_c(superkingdom, kingdom, phylum, class, order,
                     family, genus, species, sep = ";")) %>%
  pull(tax)

writeXStringSet(setNames(clustered_seqs, dada2), dada2_file)

# --- rdp
rdp = db_clus  %>%
  mutate_all(~str_replace_na(.x)) %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  mutate(tax = str_c("Root", superkingdom, kingdom, phylum, class, order,
                     family, genus, species, sep = ";")) %>%
  select(taxid, tax) %>%
  mutate(tax = str_c(taxid, tax, sep = " ")) %>%
  pull(tax)

writeXStringSet(setNames(clustered_seqs, rdp), rdp_file)

# --- mothur
mothur = db_clus  %>%
  mutate_all(~str_replace_na(.x)) %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  mutate(tax = str_c(superkingdom, kingdom, phylum, class, order,
                     family, genus, species, sep = ";")) %>%
  select(taxid, tax)

writeXStringSet(clustered_seqs, mothur_file)
write_tsv(mothur, mothur_tax, col_names = FALSE)
