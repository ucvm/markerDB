# ----------------------------------
#  Trim sequences
# ----------------------------------


library(rtracklayer)
library(Biostrings)
library(dplyr)
library(stringr)
library(purrr)

search_file = snakemake@input$hits
seq_file = snakemake@input$seqs
max_width = snakemake@config$max_width
min_width = snakemake@config$min_width
return_untrimmed = snakemake@config$return_untrimmed
outfile = snakemake@output$seqs
marker = snakemake@config$marker
search_type = snakemake@config$search_type


# for its2
model_names = list(
	"ITS2" = list(
		cmscan = list(upstream_name = "5_8S_rRNA",
								downstream_name = "LSU_rRNA_eukarya"),
		barrnap = list(upstream_name = "5_8S_rRNA",
								downstream_name = "28S_rRNA")),
	"18S" = list(
		cmscan = "SSU_rRNA_eukarya",
		barrnap = "18S_rRNA"
	)
)


load_functions = list(

	# function to load cmscan output
	cmscan = function(file) {
		
		cols = c(
			"idx", "target_name","target_accession","query_name","query_accession",
			"clan_name","mdl","mdl_from","mdl_to","seq_from","seq_to","strand","trunc",
			"pass","gc","bias","score","E_value","inc","olp","anyidx","afrct1","afrct2",
			"winidx","wfrct1","wfrct2","description_of_target"
		)
		
		cm = readr::read_table2(file, comment = "#", col_names = cols)
	
		cm = cm %>% 
			group_by(query_name, target_name) %>%
			filter(score == max(score)) %>% ungroup() %>% 
			mutate(new_seq_from = if_else(strand == "-", seq_to, seq_from),
						 new_seq_to = if_else(strand == "-", seq_from, seq_to))
		
		cm %>% 
			dplyr::select(
				seqnames = query_name,
				start = new_seq_from,
				end = new_seq_to,
				strand = strand,
				Name = target_name
			) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
	},
	
	# Function to load barrnap output
	barrnap = function(file) {
		rtracklayer::import(file, format = "gff")
	}
)


# function to load sequences
load_seqs = function(seq_file) {
	seqs = readDNAStringSet(seq_file)
	names(seqs) = str_extract(names(seqs), "^\\S+")
	return(seqs)
}

# This function takes a granges object for a single sequnece and extracts
# the start and end coordinates for trimming
get_coords = function(gr, upstream_name, downstream_name) {
		start_coord = end(gr[gr$Name == upstream_name])
		end_coord = start(gr[gr$Name == downstream_name])
		if (is_empty(start_coord)) start_coord = NA
		if (is_empty(end_coord)) end_coord = NA
		return(list(start = start_coord, end = end_coord))
}

trim_functions = list(
	# Function for processing 18S sequences
	"18S" = function(gr, seqs, names) {
		
		name = names[[1]]
		
		gr = gr[gr$Name == name]
		
		if (length(gr) == 0) {
			stop("No 18S regions detected on sequences.  You may want to check your search parameter.")
		}
		
		seqs_filtered = seqs[gr]
		names(seqs_filtered) = seqnames(gr)
		
		return(seqs_filtered)
		
	},
	
	# Trim sequences to contain only the ITS region, using a granges object
	# Returns trimmed sequences from seqs that are also in seqnames(gr)
	ITS2 = function(gr, seqs, names) {
		if (length(names) != 2) {
			stop("Please provide the names of the ITS2 flanking regions.  Should be a vector of length 2")
		}
		upstream_name = names[[1]]
		downstream_name = names[[2]]
		
		split_gr = split(gr, seqnames(gr))
		coords = lapply(split_gr, get_coords, upstream_name, downstream_name) 
	
		trimmed = 
			imap(coords, function(x, y) {
				subseq(seqs[[y]], start = x$start, end = x$end)
				}) %>% 
			DNAStringSet(use.names = TRUE)
		
		return(trimmed)
	}
)



# Run it ------------------------------------------------------------------


gr = load_functions[[search_type]](search_file)
seqs = load_seqs(seq_file)
trimmed = trim_functions[[marker]](gr, seqs, model_names[[marker]][[search_type]])

if (return_untrimmed & marker == "ITS2") {
	trimmed = c(trimmed, seqs[!names(seqs) %in% names(trimmed)])
}


# --- filter to proper lengths and keep only unqiue
message("Removing ", sum(width(trimmed) <= min_width), " sequences less than ", min_width, " bp" )
trimmed = trimmed[width(trimmed) > min_width]

message("Removing ", sum(width(trimmed) >= max_width), " sequences greater than ", max_width, " bp" )
trimmed = trimmed[width(trimmed) < max_width]

trimmed_uniq = unique(trimmed)
message("Removed ", length(trimmed) - length(trimmed_uniq), " non-unique sequences")


# Write it out ------------------------------------------------------------

writeXStringSet(trimmed_uniq, outfile)
message("Wrote ", length(trimmed_uniq), " sequences to file")



