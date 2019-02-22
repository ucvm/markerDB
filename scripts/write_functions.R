# ----------------------------------
#  Write out database formats
# ----------------------------------

# --- dada2


write_functions = list(
	dada2 = function(db, seqs, outdir, align = NULL) {
		
		outfasta = file.path(outdir,"dada2.fasta")
		
		dada2 = db %>%
			dplyr::mutate_all(~stringr::str_replace_na(.x)) %>%
			dplyr::mutate(species = stringr::str_replace_all(species, " ", "_")) %>%
			dplyr::mutate(tax = stringr::str_c(superkingdom, kingdom, phylum, class, order,
																				 family, genus, species, sep = ";")) %>%
			dplyr::pull(tax)
		
		Biostrings::writeXStringSet(setNames(seqs, dada2), outfasta)
		
		return(outfasta)
	},
	
	rdp = function(db, seqs, outdir, align = NULL) {
		
		outfasta = file.path(outdir, "rdp.fasta")
		
		rdp = db %>%
			dplyr::mutate_all(~stringr::str_replace_na(.x)) %>%
			dplyr::mutate(species = stringr::str_replace_all(species, " ", "_")) %>%
			dplyr::mutate(tax = stringr::str_c("Root", superkingdom, kingdom, phylum, class, order,
																				 family, genus, species, sep = ";")) %>%
			dplyr::select(accn, tax) %>%
			dplyr::mutate(tax = stringr::str_c(accn, tax, sep = " ")) %>%
			dplyr::pull(tax)
		
		Biostrings::writeXStringSet(setNames(seqs, rdp), outfasta)
		
		return(outfasta)
	},
	
	mothur = function(db, seqs, outdir, align = NULL) {
		
		outfasta = file.path(outdir, "mothur.fasta")
		outtaxa = file.path(outdir, "mothur.tax")
		outaln = file.path(outdir, "mothur.aln")
		
		mothur = db %>%
			dplyr::mutate_all(~stringr::str_replace_na(.x)) %>%
			dplyr::mutate(species = stringr::str_replace_all(species, " ", "_")) %>%
			dplyr::mutate(tax = stringr::str_c(superkingdom, kingdom, phylum, class, order,
																				 family, genus, species, sep = ";")) %>%
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




