
configfile: "config.yaml"

# set api key 
os.environ["ENTREZ_KEY"] = config["ncbi_api"]

org = config["organism"]
outdir = config["out_directory"]

rule all:
    input:
        "{outdir}/{org}_seqinfo.tsv".format(outdir = outdir, org = org),
        "{outdir}/{org}_raw.fasta".format(outdir = outdir, org = org),
        "{outdir}/{org}_barrnap.gff".format(outdir = outdir, org = org)


rule its2_search:
    output:
        seqinfo = "{outdir}/{org}_seqinfo.tsv".format(outdir = outdir, org = org),
        seqs = "{outdir}/{org}_raw.fasta".format(outdir = outdir, org = org)
    conda: 
        "envs/its2_search.yaml"
    script:
        "scripts/its2_search.R"

rule barrnap:
    input:
        rules.its2_search.output.seqs
    output:
        "{outdir}/{org}_barrnap.gff".format(outdir = outdir, org = org)
    conda:
        "envs/barrnap.yaml"
    threads: 
        config["threads"]
    shell:
        """
        barrnap --kingdom euk --threads {threads} {input} > {output}
        """

