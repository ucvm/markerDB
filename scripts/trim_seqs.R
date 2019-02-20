# ----------------------------------
#  Trim sequences
# ----------------------------------


library(rtracklayer)
library(Biostrings)
library(dplyr)
library(stringr)
library(purrr)
library(readr)


# input files
search_file = snakemake@input$hits
seq_infile = snakemake@input$seqs
taxa_infile = snakemake@input$taxa

# -- parameters
max_width = snakemake@config$max_width
min_width = snakemake@config$min_width
return_untrimmed = snakemake@config$return_untrimmed
marker = snakemake@config$marker

# -- output files
seqs_outfile = snakemake@output$seqs_final
seqs_outfile_nr = snakemake@output$seqs_final_nr
taxa_outfile = snakemake@output$taxa_final
taxa_outfile_nr = snakemake@output$taxa_final_nr


# - names in the hmm model 
model_names = list(
	"ITS2" = list(upstream_name = "5_8S_rRNA",
								downstream_name = "LSU_rRNA_eukarya"),
	"18S" = "SSU_rRNA_eukarya"
)


# Functions ---------------------------------------------------------------

# -- read in the cmscan results
read_cm = function(file) {
	cols = c(
		"idx", "target_name","target_accession","query_name","query_accession",
		"clan_name","mdl","mdl_from","mdl_to","seq_from","seq_to","strand","trunc",
		"pass","gc","bias","score","E_value","inc","olp","anyidx","afrct1","afrct2",
		"winidx","wfrct1","wfrct2","description_of_target"
	)
	
	cm = read_table2(file, comment = "#", col_names = cols)
	return(cm)
}


# function to load sequences
read_seqs = function(seq_file) {
	seqs = readDNAStringSet(seq_file)
	names(seqs) = str_extract(names(seqs), "^\\S+")
	return(seqs)
}

	
# --- Clean up the cmscan results
#   Cases to account for:
#    - flipped 5.8S and LSU (for ITS2)
#    - multiple hits
#    - multiple hits but on the opposite strand
#    - hits on the negative strand


clean = function(cm, seqs, names, max_hits) {

	# keep only 5.8S and LSU hits
	cm = cm %>% 
		filter(target_name %in% names)
	
	if (nrow(cm) == 0) {
		stop("No rRNA regions detected on sequences.  You may want to check your search parametesr.")
	}
	
	# keep only max score (will leave duplicate hits)
	cm = cm %>% 
		group_by(query_name, target_name) %>%
		filter(score == max(score)) %>% ungroup()
	
	# discard sequencs with more than max_hits (should be no more than 2 for ITS2 and 1 for 18S)
	discards = count(cm, query_name, sort = TRUE) %>% 
		filter(n > max_hits) %>%  pull(query_name)
	cm = filter(cm, !query_name %in% discards)

		
	# --- split by strand, if needed
	
	if ("-" %in% cm$strand) {
		by_strand = cm %>%
			split(.$strand)
	
		# identify hits from same query on opposite strands
		bad_strands = reduce(by_strand, ~intersect(.x$query_name, .y$query_name))
		# add to discard list
		discards = c(discards, bad_strands)
		
		# and filter out from hits dataframes
		by_strand = by_strand %>% map(~filter(.x, !query_name %in% bad_strands))
		
		# reverse complement the reverse strand sequenecs
		seqs[by_strand$`-`$query_name] <- reverseComplement(seqs[by_strand$`-`$query_name])
		
		# now flip the hit coords on the - hits and add these back in
		cm = by_strand$`-` %>% 
			dplyr::rename(seq_to = seq_from, seq_from = seq_to ) %>% 
			bind_rows(by_strand$`+`)
	}
	
	# sanity check: from is greater than 2
	if (any(cm$seq_from > cm$seq_to)) {
		stop("From is greater than to!  Shouldn't be here!")
	}
	
	# get rid of discards
	seqs = seqs[!names(seqs) %in% discards]
	seq_info = tibble(query_name = names(seqs), width = width(seqs))
	
	# add in sequence widths and get rid of extra columns
	cm = left_join(cm, seq_info, by = "query_name") %>% 
		dplyr::select(query_name, target_name, seq_from, seq_to, width)
	
	return(list(cm = cm, seqs = seqs))
	
}

# --- a function for each potential marker to generate trimming coordinates

coord_functions = list(
	"ITS2" = function(targets, from, to, width, name) {
		
		if (length(name) != 2) {
			stop("Please provide the names of the ITS2 flanking regions.  Should be a vector of length 2")
		}
		
		# --- Get the trimming coordinates for ITS2 region using the 
		#     upstream (5.8S) and downstream (LSU) hits
		
		upstream_name = name[[1]]
		downstream_name = name[[2]]
		
		# start at the end of the 5.8S
		start_coord = to[targets == upstream_name]
		# end at the start of the 28S (LSU)
		end_coord = from[targets == downstream_name]
		
		# If we have sequneces with one or the other subunit, it's better to 
		# explicitly set our coordinates to the start and end of the sequences
		# here so no errors later
		if (is_empty(start_coord)) start_coord = 1
		if (is_empty(end_coord)) end_coord = unique(width)
		
		# Edge case: when the end of the 5.8S is the end of the sequnce.
		#  Here we move the start coordinate back by one so the trimming
		#  will work but will result in a short sequence that gets trimmed
		if (start_coord == end_coord) start_coord = start_coord - 1
		
		
		coords = list(start = start_coord, 
									end = end_coord)
		return(coords)
	},
	
	"18S" = function(targets, from, to, width, name) {
		
		return(list(start = from, end = to))
		
	}
)

# --- Trimming function

trim = function(cm, seqs, names) {

	if (marker == "18S") {
		max_hits = 1
	} else {
		max_hits = 2
	}
	
	cleaned = clean(cm , seqs, names, max_hits)
	
	coords = cleaned$cm %>% 
		split(.$query_name) %>% 
		imap(~coord_functions[[marker]](.x$target_name, .x$seq_from, .x$seq_to, .x$width, .y)) 
	
	# and there will still be some with backwards hits, so trash these, seqs too
	backwards = names(coords[map_lgl(coords, ~.x$start > .$end | .x$start == 0 | .x$end == 0)])
	cleaned$seqs = cleaned$seqs[!names(cleaned$seqs) %in% backwards]
	coords = coords[!names(coords) %in% backwards]
	
	# now trim
	trimmed = imap(coords, function(x, y) {
		subseq(cleaned$seqs[[y]], start = x$start, end = x$end) 
	}) %>% DNAStringSet(use.names = TRUE)
	
	return(list(trimmed = trimmed, full = cleaned$seqs))

	}




# Run it ------------------------------------------------------------------


# load the files
cm = read_cm(search_file)
seqs = read_seqs(seq_infile)
taxa = read_tsv(taxa_infile)


# and now trim
trimmed = trim(cm, seqs, model_names[[marker]])

if (return_untrimmed & marker == "ITS2") {
	final_seqs = c(trimmed$trimmed, trimmed$full[!names(trimmed$full) %in% names(trimmed$trimmed)])
} else {
	final_seqs = trimmed$trimmed
}


# --- filter to proper lengths and keep only unqiue
message("Removing ", sum(width(final_seqs) <= min_width), " sequences less than ", min_width, " bp" )
final_seqs = final_seqs[width(final_seqs) > min_width]

message("Removing ", sum(width(final_seqs) >= max_width), " sequences greater than ", max_width, " bp" )
final_seqs = final_seqs[width(final_seqs) < max_width]

final_seqs_nr = unique(final_seqs)
message("Removed ", length(final_seqs) - length(final_seqs_nr), " non-unique sequences")

# and cleaned up taxonomy
final_tax = filter(taxa, accn %in% names(final_seqs))
final_tax_nr = filter(final_tax, accn %in% names(final_seqs_nr))


# Write it out ------------------------------------------------------------

writeXStringSet(final_seqs, seqs_outfile)
writeXStringSet(final_seqs_nr, seqs_outfile_nr)
write_tsv(final_tax, taxa_outfile)
write_tsv(final_tax_nr, taxa_outfile_nr)


message("Wrote ", length(final_seqs), " sequences to ", seqs_outfile)
message("Wrote ", length(final_seqs_nr), " non-redundant sequences to ", seqs_outfile_nr)



