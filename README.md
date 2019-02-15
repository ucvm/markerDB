# Build marker gene databases

## Overview

markerDB is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to build databases of a marker gene sequence for a given taxonomy.  Currently it supports ITS2 and 18S sequences.  The pipeline retrieves sequences from NCBI annotated with the provided marker, identifies the correct region (with an hmm or covariance model) and formats them for common pipelins like RDP, dada2 and mothur.

## Install :computer:

The only dependency to run the pipeline is an up-to-date conda install with a recent version of Snakemake.  If you don't already use conda then you can follow the instructions on the [Conda](https://conda.io/docs/user-guide/install/index.html) website.  Make sure to install the *Python 3* version of Miniconda.

Once conda is working then all you need to do is install Snakemake with:

```
conda install -c bioconda -c conda-forge snakemake
```

If you don't want to install Snakemake in your default environment then first create a new environment with Python 3 and Snakemake installed.

```
conda create -n markerDB python=3 snakemake
```

Then activate it:

```
source activate markerDB
```

Then you'll want to make sure you have an [NCBI account](https://www.ncbi.nlm.nih.gov/account/) from which you'll get an API Key on the account setting page.

## Running the pipeline :gear: 

First clone the repository to a location desired:

```
git clone https://github.com/ucvm/markerDB.git
```

### Configuration

Edit the `config.yaml` file with a file editor of your choice and adjust the following parameters as you see fit.

* *organism*: This is added to the NCBI search string to restrict the search to a taxonomy of interest.

* *marker*: The marker to search for.  One of "ITS2" or "18S"

* *search_type*: The search program to use.  cmscan is from [infernal](http://eddylab.org/infernal/) and uses covariance models that include secondary structure.  It's very accurate but slower.  [barrnap](https://github.com/tseemann/barrnap) uses nhmmer from HMMER 3 and is much faster but may miss the endpoints by a few basepairs (although this is likely to matter little for taxonomy assignment).

* *out_directory*: Output directory to store the files.  Will be created if it doesn't exist

* *max_width*: The final sequences will be trimmed to be no more than this many basepairs long.

* *min_width*: As above but minimum width

* *return_trimmed*: For ITS2 only.  Should sequences that didn't have a 5.8S or 28S be returned.  Default is TRUE.  See below for details

* *ncbi_api*: Because this pipeline searchs the NCBI nucleotide database and downloads thousands of sequences it's required to have an account and API key (see above on how to get one).

* *threads*: The number of CPU threads that Snakemake can use to run jobs.

### Run it

Now run it with 2 threads:

```
snakemake -j 2 --use-conda
```

The output databases will be named by their format (see Details below) with `.fasta`.  The other files can be deleted to save space but usually aren't too big and can save time if you want to rerun some of the steps without re-downloading all the sequences again.

The first run will take a little bit longer as Snakemake will build the conda environments.  

### Compute requirements

Mostly time :hourglass:.  `cmscan` is the main time consuming process and can be sped up to some degree by providing more threads.  `barrnap` is a fair bit quicker.  The pipeline doesn't require excesive memory and can likely be run on a decently powered desktop in a few hours. But this depends strongly on the number of sequences being downloaded. Snakemake has excellent support for compute clusters so if you have access to one this is recommended.  Please see the Snakemake documentation for more details on how to do this.  Also note that this will run on MacOS or Linux only as with the majority of bioinformatics applications.


## Details :mag:

### Step 1: Get sequences

This step utilizes the `rentrez` R package to search NCBI for sequences annoated as ITS2 for your given organism.  The exact search terms are:

```
ITS2: {organism}[ORGN] (ITS2 OR internal transcribed spacer 2)
18S:  {organsim}[ORGN] (18S OR SSU rRNA OR SSU OR rRNA)
```

where `{organsism}` is replaced with whatever is specified in the config file.

Taxonomy for each of the sequences is retrieved from NBCI Taxonomy with the `taxize` R package and then the fasta sequences are downloaded.

### Step 2: Trim sequences

For the 18S the sequences are trimmed to the region identified in the search.  However, the ITS2 is more complicated as no good sequence models exist that capture a wide range of diversity due to the divergence of this region.  Relying soley on the annotion in NCBI is also not feasable as many of the ITS2 sequences will also contain the other ribosomal genes, especially the 5.8S and partial 28S genes.  

After the additional rRNA genes are identified the sequences are trimmed to this region.  If no additional rRNA genes were found the sequence is returned untouched if the `return_untrimmed` parameter is set.  Otherwise these are discarded.  It is suggested to leave this set to `TRUE` as many ITS2 sequences have been deposited to NCBI already trimmed to that region.

The sequences are then clipped to at most `max_width` and at minimum `min_width`.  Also any redundant sequences that may be present so the database is unique.

### Step 4: Alignment

Although not required for many taxonomic classifiers, like RDP, some pipelines do require an alignment the database is aligned using mafft.

### Step 5:  Format database

Currently markerDB outputs 3 common formats used for assigning taxonomy.  

* dada2: dada2's `assignTaxonomy` function
* RDP: to train a custom RDP database with the `rRDP` Bioconductor package (this is pretty much the same as `assignTaxonomy`)
* mothur: a fasta file and paired mothur taxonomy file.  Works with the [Nemabiome](https://www.nemabiome.ca/) pipeline.  An alignment is also written out that should work with mothur, but can also be used for other pipelines as needed.

The seqinfo.tsv file is a tab-separated text file that contains the taxonomy information for all sequences and can be easily imported into R or Python or (if you must) be opened in Excel.

Files appropriate for [IDTAXA](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0521-5) will be coming shortly.


## TODO :hammer:

* Add update step to avoid redownloading existing sequences
  - this will likely entail an ids file that lists all searched sequences
* Some sort of filtering on per species basis.  There are some clear outliers that need to be cleaned out.

## For the pros :trophy:

For doing it your own way the best bet is to fork the repository and edit away.  This way your changes will be preserved so you can publishe them properly!  Pull requests are welcome but need to be broadly useful and well tested.

Other typos/errors can also be corrected with pull requests or by submitting an issue.







