
configfile: "config.yaml"

# set api key 
os.environ["ENTREZ_KEY"] = config["ncbi_api"]

org = config["organism"]
outdir = config["out_directory"]

rule all:
    input:
        "{outdir}/{org}_seqinfo.tsv".format(outdir = outdir, org = org),
        "{outdir}/{org}_raw.fasta".format(outdir = outdir, org = org),
        "{outdir}/{org}_cmscan.out".format(outdir = outdir, org = org),
        "{outdir}/{org}_barrnap.out".format(outdir = outdir, org = org),
        "{outdir}/{org}_trimmed.fasta".format(outdir = outdir, org = org),
        "{outdir}/{org}_rdp.fasta".format(outdir = outdir, org = org),
        "{outdir}/{org}_dada2.fasta".format(outdir = outdir, org = org),
        "{outdir}/{org}_mothur.fasta".format(outdir = outdir, org = org),
        "{outdir}/{org}_mothur.tax".format(outdir = outdir, org = org),
        "{outdir}/{org}_mothur.aln".format(outdir = outdir, org = org)


rule sequence_search:
    output:
        seqinfo = "{outdir}/{org}_seqinfo.tsv".format(outdir = outdir, org = org),
        seqs = "{outdir}/{org}_raw.fasta".format(outdir = outdir, org = org)
    conda: 
        "envs/its2_search.yaml"
    script:
        "scripts/its2_search.R"

rule cmscan:
    input: 
        rules.sequence_search.output.seqs
    output:
        "{outdir}/{org}_cmscan.out".format(outdir = outdir, org = org)
    params: 
        models = "cm_model/marker_cm_database.cm"
    conda:
        "envs/cmscan.yaml"
    threads: 
        config["threads"]
    shell:
        """
        cmscan --rfam --cut_ga --nohmmonly --tblout {output} --fmt 2 --cpu {threads} {params.models} {input} > /dev/null 2>&1
        """
        
rule barrnap:
	input:
		rules.sequence_search.output.seqs
	output:
		"{outdir}/{org}_barrnap.out".format(outdir = outdir, org = org)
	conda:
		"envs/barrnap.yaml"
	threads: 
		config["threads"]
	shell:
		"""
		barrnap -k euk -t 8 -l 0.25 -r 0.1 --quiet {input} > {output}
		"""

rule trim_seqs:
    input:
        #hits = rules.cmscan.output,
        hits = rules.barrnap.output,
        seqs = rules.sequence_search.output.seqs
    output:
        seqs = "{outdir}/{org}_trimmed.fasta".format(outdir = outdir, org = org)
    conda: 
        "envs/trim_seqs.yaml"
    script:
        "scripts/trim_seqs.R"

rule write_seqs:
    input:
        seqs = rules.trim_seqs.output.seqs,
        info = rules.sequence_search.output.seqinfo
    output:
        rdp = "{outdir}/{org}_rdp.fasta".format(outdir = outdir, org = org),
        dada2 = "{outdir}/{org}_dada2.fasta".format(outdir = outdir, org = org),
        mothur = "{outdir}/{org}_mothur.fasta".format(outdir = outdir, org = org),
        mothur_tax = "{outdir}/{org}_mothur.tax".format(outdir = outdir, org = org)
    conda:
        "envs/write_seqs.yaml"
    script:
        "scripts/write_seqs.R"
        
rule cluster_by_species:
    input:
        seqs = rules.trim_seqs.output.seqs,
        info = rules.sequence_search.output.seqinfo
    output:
        rdp = "{outdir}/{org}_rdp_clustered.fasta".format(outdir = outdir, org = org),
        dada2 = "{outdir}/{org}_dada2_clustered.fasta".format(outdir = outdir, org = org),
        mothur = "{outdir}/{org}_mothur_clustered.fasta".format(outdir = outdir, org = org),
        mothur_tax = "{outdir}/{org}_mothur_clustered.tax".format(outdir = outdir, org = org)
    conda:
        "envs/cluster.yaml"
    log:
        "logs/cluster_by_species.log"
    script:
        "scripts/cluster_by_species.R"

rule align_full:
    input:
        rules.write_seqs.output.mothur
    output:
        "{outdir}/{org}_mothur.aln".format(outdir = outdir, org = org)
    conda:
        "envs/align.yaml"
    threads:
        config["threads"]
    log:
        "logs/align_full.log"
    shell:
        """
        mafft --auto --reorder --thread {threads} {input} > {output} 2> {log}
        """

rule align_clustered:
    input:
        rules.cluster_by_species.output.mothur
    output:
        "{outdir}/{org}_mothur_clustered.aln".format(outdir = outdir, org = org)
    conda:
        "envs/align.yaml"
    threads:
        config["threads"]
    log:
        "logs/align_clustered.log"
    shell:
        """
        mafft --auto --reorder --thread {threads} {input} > {output} 2> {log}
        """


