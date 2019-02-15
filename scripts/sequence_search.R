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
marker = snakemake@config$marker



# make sure we've got a supported marker
if (!marker %in% c("ITS2", "18S")) {
	stop(marker, " not supported")
}

search_text = list(
	"ITS2" = "(ITS2 OR internal transcribed spacer 2)",
	"18S" = "18S"
)

term = glue("{organism}[ORGN] {search_text[[marker]]}")

# search nucleotide using provided organism -------------------------------

run_search = function(term) {
	message("Search NCBI with the following search term ", term)
	search = entrez_search(db = "nucleotide", term = term, use_history = TRUE)
	message("Found ", search$count, " hits")
	return(search)
}

# Get summaries -----------------------------------------------------------
# needed so we can add taxonomy
get_summaries = function(search) {
	message("Getting summaries...")
	summaries = seq(0, search$count, by = 200) %>%
	  map(~entrez_summary(db = "nucleotide",
	                      web_history = search$web_history,
	                      retmax = 200,
	                      retstart = .x)) %>%
	  flatten()
}

# turn search summaries into a dataframe
summaries_to_dataframe = function(summaries) {
	map_dfr(summaries, 
					~tibble(
						uid = as.character(.x[["uid"]]),
						taxid = as.character(.x[["taxid"]])))
}

# add accension numbers
add_acc = function(sum_df) {
	ids = split(sum_df$uid, ceiling(seq_along(sum_df$uid) / 200))
	accns = map(ids, ~entrez_fetch("nucleotide", .x, rettype = "acc"))
	accns = str_c(accns, collapse = "\n")
	accns = str_split(accns, "\n")[[1]]
	accns = discard(accns, ~.x == "")
	
	if (length(accns) != nrow(sum_df)) {
	  stop("Number of accession ids retrieved does not match the number of search hits")
	}
	
	return(mutate(sum_df, accn = accns))
}


# Get taxonomic classifications -------------------------------------------

add_taxonomy = function(summaries, sum_df) {
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
	
	
	# filter out bad taxonomies, species without proper identification,
	# contains sp. and numbers.  Also filters out taxa that are NA
	sum_df = sum_df %>% 
		filter(!str_detect(sum_df$species, "sp\\.|\\d+"),
					 !is.na(genus),
					 !is.na(family),
					 !is.na(order),
					 !is.na(class),
					 !is.na(phylum),
					 !is.na(kingdom),
					 !is.na(superkingdom))
	
	
	return(sum_df)

}


# Get fasta sequences -----------------------------------------------------

retrieve_and_write_sequences = function(ids, outfile) {
	message("Retrieving sequences...")
	ids = split(ids, ceiling(seq_along(ids) / 50))
 	seqs = map(ids, ~entrez_fetch("nucleotide", 
 																id = .x, 
 																rettype = "fasta", 
 																retmode = "text")) %>% 
 		flatten()
	
 	message("Writing fasta of hits")
 	if (file.exists(outfile)) {
 		warning(outfile, " exists, removing..")
 		file.remove(outfile)
 	}
 	
	seqs %>%
  	walk(~write_file(.x, outfile, append = TRUE))
}



# Run it ------------------------------------------------------------------

search = run_search(term)
summaries = get_summaries(search)
sum_df = summaries_to_dataframe(summaries)
sum_df = add_acc(sum_df)
sum_df = add_taxonomy(summaries, sum_df)
retrieve_and_write_sequences(sum_df$uid, snakemake@output$seqs)

message("Writing taxonomy info")
write_tsv(sum_df, snakemake@output$seqinfo)



