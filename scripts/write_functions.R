# ----------------------------------
#  Write out database formats
# ----------------------------------

# --- dada2


get_ranks = function(db) {
	all_ranks = c("superkingdom", "kingdom", "phylum", "class", 
								"order", "family", "genus", "species")
		
	ranks = colnames(db)[colnames(db) %in% all_ranks]
	return(ranks)
}

clean_db = function(db) {
	db = dplyr::mutate_all(db, ~stringr::str_replace_na(.x)) 
	if ("species" %in% colnames(db)) {
		db = dplyr::mutate(db, species = stringr::str_replace_all(species, " ", "_"))
	}	
	return(db)
}


write_functions = list(
	dada2 = function(db, seqs, outdir, align = NULL) {
		
		outfasta = file.path(outdir,"dada2.fasta")
		ranks = get_ranks(db)
		
		dada2 = db %>% 
			clean_db() %>% 
			dplyr::select(ranks) %>% 
			tidyr::unite("tax", sep = ";") %>% 
			dplyr::pull(tax)
	
		Biostrings::writeXStringSet(setNames(seqs, dada2), outfasta)
		
		return(outfasta)
	},
	
	rdp = function(db, seqs, outdir, align = NULL) {
		
		outfasta = file.path(outdir, "rdp.fasta")
		ranks = get_ranks(db)
		
		rdp = db %>% 
			clean_db() %>% 
			mutate(root = "Root") %>% 
			dplyr::select(accn, root, ranks) %>% 
			tidyr::unite("tax", -accn, sep = ";") %>% 
			dplyr::mutate(tax = stringr::str_c(accn, tax, sep = " ")) %>%
			dplyr::pull(tax)
		
		Biostrings::writeXStringSet(setNames(seqs, rdp), outfasta)
		
		return(outfasta)
	},
	
	mothur = function(db, seqs, outdir, align = NULL) {
		
		outfasta = file.path(outdir, "mothur.fasta")
		outtaxa = file.path(outdir, "mothur.tax")
		outaln = file.path(outdir, "mothur.aln")
		
		ranks = get_ranks(db)
		
		mothur = db %>%
			clean_db() %>% 
			dplyr::select(accn, ranks) %>% 
			tidyr::unite("tax", -accn, sep = ";") %>% 
			dplyr::mutate(tax = stringr::str_c(tax, ";")) %>%
			dplyr::select(accn, tax)
		
		Biostrings::writeXStringSet(seqs, outfasta)
		readr::write_tsv(mothur, outtaxa, col_names = FALSE)
		file_list = list(outfasta, outtaxa)
		if (!is.null(align)) {
			Biostrings::writeXStringSet(align, outaln)
			file_list = c(outaln, file_list)
		}
		
		return(file_list)
	}
)




