
configfile: "config.yaml"

# set api key 
os.environ["ENTREZ_KEY"] = config["ncbi_api"]

#org = config["organism"]
outdir = config["out_directory"]

rule all:
    input:
        "{outdir}/db/seqs.fasta".format(outdir = outdir),
        "{outdir}/db/seqs_nr.fasta".format(outdir = outdir),
        "{outdir}/db/taxonomy.txt".format(outdir = outdir),
        "{outdir}/db/taxonomy_nr.txt".format(outdir = outdir)


rule sequence_search:
    output:
        taxa = "{outdir}/raw/raw_taxa.tsv".format(outdir = outdir),
        seqs = "{outdir}/raw/raw_seqs.fasta".format(outdir = outdir)
    conda: 
        "envs/sequence_search.yaml"
    script:
        "scripts/sequence_search.R"

rule cmscan:
    input: 
        rules.sequence_search.output.seqs
    output:
        "{outdir}/raw/cmscan.out".format(outdir = outdir)
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
        
rule trim_seqs:
    input:
        hits = rules.cmscan.output,
        seqs = rules.sequence_search.output.seqs,
        taxa = rules.sequence_search.output.taxa
    output:
        seqs_final = "{outdir}/db/seqs.fasta".format(outdir = outdir),
        seqs_final_nr = "{outdir}/db/seqs_nr.fasta".format(outdir = outdir),
        taxa_final = "{outdir}/db/taxonomy.txt".format(outdir = outdir),
        taxa_final_nr = "{outdir}/db/taxonomy_nr.txt".format(outdir = outdir)
    conda: 
        "envs/trim_seqs.yaml"
    script:
        "scripts/trim_seqs.R"

rule write_seqs:
    input:
        seqs = rules.trim_seqs.output.seqs_final,
        taxa = rules.trim_seqs.output.taxa_final,
    output:
        rdp = "{outdir}/formats/rdp.fasta".format(outdir = outdir),
        dada2 = "{outdir}/formats/dada2.fasta".format(outdir = outdir),
        mothur = "{outdir}/formats/mothur.fasta".format(outdir = outdir),
        mothur_tax = "{outdir}/formats/mothur.tax".format(outdir = outdir)
    conda:
        "envs/write_seqs.yaml"
    script:
        "scripts/write_seqs.R"
        
        
rule write_seqs_nr:
    input:
        seqs = rules.trim_seqs.output.seqs_final_nr,
        taxa = rules.trim_seqs.output.taxa_final_nr
    output:
        rdp = "{outdir}/formats/rdp_nr.fasta".format(outdir = outdir),
        dada2 = "{outdir}/formats/dada2_nr.fasta".format(outdir = outdir),
        mothur = "{outdir}/formats/mothur_nr.fasta".format(outdir = outdir),
        mothur_tax = "{outdir}/formats/mothur_nr.tax".format(outdir = outdir)
    conda:
        "envs/write_seqs.yaml"
    script:
        "scripts/write_seqs.R"
        
        
# rule cluster_by_species:
#     input:
#         seqs = rules.trim_seqs.output.seqs,
#         info = rules.sequence_search.output.seqinfo
#     output:
#         rdp = "{outdir}/{org}_rdp_clustered.fasta".format(outdir = outdir, org = org),
#         dada2 = "{outdir}/{org}_dada2_clustered.fasta".format(outdir = outdir, org = org),
#         mothur = "{outdir}/{org}_mothur_clustered.fasta".format(outdir = outdir, org = org),
#         mothur_tax = "{outdir}/{org}_mothur_clustered.tax".format(outdir = outdir, org = org)
#     conda:
#         "envs/cluster.yaml"
#     log:
#         "logs/cluster_by_species.log"
#     script:
#         "scripts/cluster_by_species.R"


rule align_full:
    input:
        rules.write_seqs.output.mothur
    output:
        "{outdir}/formats/mothur.aln".format(outdir = outdir)
    conda:
        "envs/align.yaml"
    threads:
        config["threads"]
    shell:
        """
        mafft --auto --reorder --thread {threads} {input} > {output}
        """



