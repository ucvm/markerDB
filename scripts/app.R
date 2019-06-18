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
library(purrr)
library(tibble)
library(glue)

Sys.setenv(R_ZIPCMD = "/usr/bin/zip")

source("write_functions.R")
db_dir = "../../Nematode_ITS2_1.0.0/"


#db_dir = snakemake@params$db_dir
#source("scripts/write_functions.R")


load_database = function(db_directory) {

	seqs = readDNAStringSet(file.path(db_directory, "db", "seqs.fasta"))
	seqs_nr = readDNAStringSet(file.path(db_directory, "db", "seqs_nr.fasta"))
	align_nr = readDNAStringSet(file.path(db_directory, "db", "seqs_nr.aln"))
	taxa = read_tsv(file.path(db_directory, "db", "taxonomy.txt"),
					col_types = "cccffffffff")
	taxa_nr = read_tsv(file.path(db_directory, "db", "taxonomy_nr.txt"),
					   col_types = "cccffffffff")
	
	
	db = tibble(accn = names(seqs), sequence = as.character(seqs)) %>% 
		left_join(taxa) %>% dplyr::select(-sequence, everything())
	db_nr = tibble(accn = names(seqs_nr), sequence = as.character(seqs_nr)) %>% 
		left_join(taxa_nr) %>% dplyr::select(-sequence, everything())
	
	return(list(db = db, nr = db_nr, seqs = seqs, seqs_nr = seqs_nr, align_nr = align_nr))
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
db_info = read_tsv(file.path(db_dir, "db", "db_info.txt"), comment = "#", 
									 col_names = c("param", "value")) %>% deframe()

db_info[3:4] = prettyNum(db_info[3:4], big.mark = ",")

db_html_list = imap(db_info, 
										~p(HTML(glue::glue('<strong>{.y}</strong>: {.x}'))))[-1]



ui = fluidPage(theme = shinytheme("yeti"),

	useShinyjs(),
	

    verticalLayout(
    	
    	titlePanel(db_info[["Title"]]),
    
    	fluidRow(
    		column(3,     		
				 h3("Summary"),
				 wellPanel(div(db_html_list),
				 					tags$small("This database was created with markerDB",
				 										 a(icon("github"), href = "https://github.com/ucvm/markerDB", target = "_blank"))
    			)
				 ),
				column(9,
					h3("Instructions"),
					p("Please use the options and table below to filter the database as required.", 
						"Once it's been filtered to your liking you can download the filtered database.",
						"Please keep track of your filtering steps and note these if you publish."),
					h4("Download formats"),
					tags$ul(
						tags$li(
							strong("dada2:"), 
							"A fasta file with taxonomy in the headers, suitable for dada2::assignTaxonomy"
							),
						tags$li(
							strong("mothur:"), 
							"A fasta file with sequence ids as headers, a taxonomy text file, and an",
							"alignment generated with mafft"
							),
						tags$li(
							strong("rdp:"), 
							"A fasta file with sequence ids and taxonomy in the headers.", 
							"Used to train a custom RDP database with the ", code("rRDP"), 
							"Biodondcutor package (which uses the same algorithm as ", code("assignTaxonomy"), "."
						),
						tags$li(
							strong("idtaxa:"), 
							"A fasta file with taxonomy in the headers, and a taxid file suitable for", 
							"IDTAXA with DECIPHER::LearnTaxa"
						),
						tags$li(
							strong("excel: "),
							"An Excel file with the same columns in the table."
						)
					)
    		)
    	),

    	hr(),
    	fluidRow(
    		column(3,
	
	     		h3("Current Selection:"),
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
									   c("dada2", "mothur", "rdp", "idtaxa", "excel"), inline = FALSE),
					em("Please select a format to download", id = "no_format_text"),
					hidden(em("The accn column must be selected for database download", id = "no_accn_text")),
					hidden(downloadButton("download", "Download"))
    			),
					hidden(
						div("You've filtered the database.  If you use the filtered version in an analyis",
								"or publication, please document what you've filtered out to maintain reproducibility",
								class = "alert alert-warning", role = "alert", id = "filter_alert")
					)
    		)
    	),
    	
    	hr(),
			h3("Database"),
			br(),
    	DT::DTOutput("database")
    )
)

server = function(input, output) {
	
	# --- choose version to show
	filtered_data = reactive({
		if (input$show_full) {
			display = database$db
		} else {
			display = database$nr
		}
		display[ , input$show_columns]
	})
	
	# --- control when to show download button
	observe({
		if (is.null(input$download_formats)) {
			shinyjs::hide("download")
			shinyjs::show("no_format_text")
			shinyjs::hide("no_accn_text")
		} else if (!"accn" %in% colnames(filtered_data())) {
			shinyjs::hide("download")
			shinyjs::hide("no_format_text")
			shinyjs::show("no_accn_text")
		} else {
			shinyjs::show("download")
			shinyjs::hide("no_format_text")
			shinyjs::hide("no_accn_text")
		}
	})
	
	# --- summary table
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
	

	output$database = DT::renderDataTable({
		
		DT::datatable(filtered_data(), style = "bootstrap", 
									class = "table-hover table-condensed table-striped table-responsive",
									filter = "top")
	})

		
	# --- filter alert
	observe({	
		if (input$show_full) {
			full = database$db
		} else {
			full = database$nr
		}
		
		if (length(input$database_rows_all) < nrow(full)) {
			shinyjs::show("filter_alert")
		} else {
			shinyjs::hide("filter_alert")
		}
	})
	
	
	# --- download
	output$download = downloadHandler(
		filename = "database.zip",
		content = function(file) {
			tmpdir = tempdir()
			db = filtered_data()[input$database_rows_all, ]
			seqs = database$seqs[db$accn]
			if (!input$show_full) {
				aln  = database$align_nr[db$accn]
			} else {
				aln = NULL
			}
			setwd(tmpdir)
			files = c()
			write_formats = names(write_functions) %in% input$download_formats
			if (!is.null(write_formats)) {
				db_files = invoke_map(write_functions[write_formats],
								   list(list(db = db, seqs = seqs, outdir = ".", align = aln)))
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
	
	
}

# Run the application 
shinyApp(ui = ui, server = server)
