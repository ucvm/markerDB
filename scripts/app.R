# ----------------------------------
#  A shiny web app to view the database
# ----------------------------------

library(Biostrings)
library(shiny)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(DT)
library(shinyjs)
library(shinythemes)

source("scripts/write_functions.R")
Sys.setenv(R_ZIPCMD = "/usr/bin/zip")

# db_dir = "../Nematode_ITS2/"
# db_title = "Nematode ITS2"

db_dir = snakemake@params$db_dir
db_title = snakemake@params$db_title

load_database = function(db_directory) {

	seqs = readDNAStringSet(file.path(db_directory, "db", "seqs.fasta"))
	seqs_nr = readDNAStringSet(file.path(db_directory, "db", "seqs_nr.fasta"))
	taxa = read_tsv(file.path(db_directory, "db", "taxonomy.txt"),
					col_types = "cccffffffff")
	taxa_nr = read_tsv(file.path(db_directory, "db", "taxonomy_nr.txt"),
					   col_types = "cccffffffff")
	
	db = tibble(accn = names(seqs), sequence = as.character(seqs)) %>% 
		left_join(taxa) %>% dplyr::select(-sequence, everything())
	db_nr = tibble(accn = names(seqs_nr), sequence = as.character(seqs_nr)) %>% 
		left_join(taxa_nr) %>% dplyr::select(-sequence, everything())
	
	return(list(db = db, nr = db_nr, seqs = seqs, seqs_nr = seqs_nr))
}

summarize_database = function(db) {
	db %>% 
		select(-accn, -uid, -taxid, -sequence) %>% 
		gather(Rank, taxa) %>% 
		group_by(Rank) %>% 
		summarise(`Number of taxa` = n_distinct(taxa)) %>% 
		mutate(Rank = str_to_title(Rank)) %>%
		arrange(factor(Rank, levels = c("Superkingdom", "Kingdom", "Phylum", "Class", "Order",
										"Family", "Genus", "Species")))
		
}

database = load_database(db_dir)

ui = fluidPage(theme = shinytheme("yeti"),

	useShinyjs(),
	

    verticalLayout(
    	
    	titlePanel(db_title),
    
    	
    	fluidRow(
    		column(3,
	     		h3("Summary"),
	     		wellPanel(
        			tableOutput("summary_tbl")
	     		)
    		),
    		
    		column(4,
    			h3("Viewing options"),
    			wellPanel(
					checkboxGroupInput("show_columns", "Select columns to display", colnames(database$nr),
							   selected = colnames(database$nr)[!colnames(database$nr) == "sequence"],
							   inline = TRUE),
					hr(),
        			checkboxInput("show_full", "Include redundant (identical) sequences", value = FALSE),
					helpText("This is not recommended unless you have a specific reason for",
							 "doing it.  The non-redundant version is prefered for taxonomy assignment")
    			)
    			   
    		),
    		
    		column(3, 
    			h3("Download"),
    			wellPanel(
					checkboxGroupInput("download_formats", "Select formats to download",
									   c("dada2", "mothur", "rdp", "excel"), inline = TRUE),
					em("Please select a format to download", id = "no_format_text"),
					hidden(em("The accn column must be selected for database download", id = "no_accn_text")),
					hidden(downloadButton("download", "Download"))
    			)
    		)
        ),
    	hr(),
        DT::DTOutput("database")
    )
)

server = function(input, output) {
	
	filtered_data = reactive({
		# --- choose version to show
		if (input$show_full) {
			display = database$db
		} else {
			display = database$nr
		}
		display[ , input$show_columns]
	})
	
	observe({
		if (is.null(input$download_formats)) {
			hide("download")
			show("no_format_text")
			hide("no_accn_text")
		} else if (!"accn" %in% colnames(filtered_data())) {
			hide("download")
			hide("no_format_text")
			show("no_accn_text")
		} else {
			show("download")
			hide("no_format_text")
			hide("no_accn_text")
		}
	})
	
	output$summary_tbl = renderTable({
		ranks = c("superkingdom", "kingdom", "phylum", "class", 
				  "order","family", "genus", "species")
		filtered_data()[input$database_rows_all, ] %>% 
			select(one_of(ranks)) %>% 
			gather(Rank, taxa) %>% 
			group_by(Rank) %>% 
			summarise(`Number of taxa` = n_distinct(taxa)) %>% 
			mutate(Rank = str_to_title(Rank)) %>%
			arrange(factor(Rank, levels = c("Superkingdom", "Kingdom", "Phylum", "Class", "Order",
											"Family", "Genus", "Species")))
	}, bordered = TRUE, spacing = "xs")
	
	output$download = downloadHandler(
		filename = "database.zip",
		content = function(file) {
			tmpdir = tempdir()
			db = filtered_data()[input$database_rows_all, ]
			seqs = database$seqs[db$accn]
			setwd(tmpdir)
			files = c()
			write_formats = names(write_functions) %in% input$download_formats
			if (!is.null(write_formats)) {
				db_files = invoke_map(write_functions[write_formats],
								   list(list(db = db, seqs = seqs, outdir = ".")))
				db_files = flatten(db_files) %>% flatten_chr()
				files = c(files, db_files)
			}
			if ("excel" %in% input$download_formats) {
				writexl::write_xlsx(db, "database.xlsx")
				files = c(files, "database.xlsx")
			}
			
			zip(file, files)

		},
		contentType = "application/zip"
	)
	
	output$database = DT::renderDataTable({

		DT::datatable(filtered_data(), style = "bootstrap", 
					  class = "table-hover table-condensed table-striped table-responsive",
					  filter = "top")
	})
		
}

# Run the application 
shinyApp(ui = ui, server = server)
