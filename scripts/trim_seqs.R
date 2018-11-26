# ----------------------------------
#  Trim sequences
# ----------------------------------


library(rtracklayer)
library(Biostrings)
library(dplyr)
library(stringr)
library(purrr)

cmscan_file = snakemake@input$cm
seq_file = snakemake@input$seqs
max_width = snakemake@config$max_width
min_width = snakemake@config$min_width
return_untrimmed = snakemake@config$return_untrimmed
outfile = snakemake@output$seqs


# Functions ---------------------------------------------------------------

extract_ITS = function(hits, seq, extend = 0) {

    start = hits$seq_to[which(hits$target_name == "5_8S_rRNA" )]
    start = start + 1

    end = hits$seq_from[which(hits$target_name == "LSU_rRNA_eukarya" )]
    end = end - 1

    if (length(end) == 0) {
      end = length(seq)
    } else {
      end = end + extend
      if (end > length(seq)) end = length(seq)
    }

    if (length(start) == 0) {
      start = 1
    } else {
      start = start - extend
      if (sign(start) == -1) start = 1
    }

    # 5.8S is at the end - no ITS - return empty DNAString
    if (start > length(seq)) {
      return(Biostrings::DNAString())
    }

    seq[start:end]
}

process_single = function(seq, id, cmscan_df, return_untrimmed) {
  #id = str_extract(id, "([^\\s]+)")
  if (id %in% cmscan_df$query_name) {
    hits = filter(cmscan_df, query_name == id)
    new_seq = extract_ITS(hits, seq)
  } else if (return_untrimmed) {
    new_seq = seq
  } else {
    new_seq = DNAString()
  }
}



# load seqs and cmscan output ---------------------------------------------


cols = c(
  "idx", "target_name","target_accession","query_name","query_accession",
  "clan_name","mdl","mdl_from","mdl_to","seq_from","seq_to","strand","trunc",
  "pass","gc","bias","score","E_value","inc","olp","anyidx","afrct1","afrct2",
  "winidx","wfrct1","wfrct2","description_of_target"
)

cm = read.table(cmscan_file, header = FALSE, col.names = cols, comment.char = "#",
                 stringsAsFactors = FALSE) %>% as_tibble() %>% 
  group_by(query_name, target_name) %>%
  filter(score == max(score)) %>% ungroup()

seqs = readDNAStringSet(seq_file)
names(seqs) = str_extract(names(seqs), "^\\S+")


# Trim --------------------------------------------------------------------


# --- trim to ITS region
trimmed = imap(as.list(seqs), ~process_single(.x, .y, cm, return_untrimmed))

# --- convert to DNAStringSet
trimmed = DNAStringSet(trimmed)

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



