#snakemake --dag | dot -Tpdf > dag.pdf
#snakemake --detailed-summary > provenance.tsv
#snakemake --cores 6 --resources mem=12
configfile: "config.yaml"

n_sim = config['n_sim']
cpu_type = config['cpu_type']
thrs = config['threads']

hg = config['hg']
bwa_indexes = [hg+".bwt", hg+".pac", hg+".amb", hg+".ann", hg+".sa"]

datadir = 'data/'
resultdir = 'results/'

import os
samples = [filename for filename in os.listdir('./'+datadir) if filename.endswith('.fastq')]
samples = set("_".join(filename.split('_')[:-1]) for filename in samples)


rule all:
    input:
        expand(resultdir+"{sample}_dedup.bai", sample=samples),
        hg.replace('fasta', 'dict'),
        expand(resultdir+"{sample}.interval", sample=samples),
    benchmark:
        "benchmarks/benchmark_rule_all_subset_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs)
    run:
        pass


rule mapping:
    """
    this maps the samples to the reference genome.
    It does append a header necessary for the GATK analysis.
    """
    input:
        bwa_indexes,
        sample1 = datadir+"{sample}_1.fastq",
        sample2 = datadir+"{sample}_2.fastq",
    output:
        outfile = resultdir+"{sample}.sam",
    params:
        reference = hg,
        library=config['library'],
        platform=config['platform'],
        name = "{sample}",
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_mapping_{sample}_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, sample=samples)
    threads: 2
    resources: mem=6
    version: 0.1
    message: "bwa is aligning the sample '{params.name}' with the reference genome"
    shell:
        "bwa mem -M -t {threads} -R \'@RG\\tID:{params.name}\\tSM:{params.name}\\tPL:{params.platform}\\tLB:{params.library}\\tPU:PE\' {params.reference} {input.sample1} {input.sample2} > {output.outfile}"


rule sort_picard:
    input:
        r = resultdir+"{sample}.sam",
    output:
        resultdir+"{sample}_sorted.bam",
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_sort_picard_{sample}_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, sample = samples)
    shell:
        "picard SortSam INPUT={input.r} OUTPUT={output}.tmp SORT_ORDER=coordinate"
        " && [ -s {output}.tmp ] && mv {output}.tmp {output}"


rule mark_duplicates:
    input:
        r=resultdir+"{sample}_sorted.bam",
    output:
        resultdir+"{sample}_dedup.bam",
    params:
        outdir=resultdir,
        metricsfile=resultdir+"{sample}_dedup.metrics.txt",
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_mark_duplicates_{sample}_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, sample = samples)
    shell:
        "picard MarkDuplicates"
        " INPUT={input.r} OUTPUT={output}.tmp METRICS_FILE={params.metricsfile}"
        " MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000"
        " && [ -s {output}.tmp ] && mv {output}.tmp {output}"


rule build_bam_index:
    input:
        r=resultdir+"{sample}_dedup.bam",
    output:
        resultdir+"{sample}_dedup.bai",
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_build_bam_{sample}_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, sample = samples)
    shell:
        "picard BuildBamIndex INPUT={input.r} OUTPUT={output}"


rule realigner_target_creator:
    input:
        seq=resultdir+"{sample}_dedup.bam",
        idx=resultdir+"{sample}_dedup.bai",
        ref=hg+'.fai',
    output:
        resultdir+"{sample}.interval",
    params:
        gatk = config['gatk'],
        #gatk='programs/gatk/GenomeAnalysisTK.jar',
        realref=hg,
    benchmark:
        "benchmarks/benchmark_realigner_{sample}_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, sample = samples)
    shell:
        "java -jar {params.gatk} -T RealignerTargetCreator -R {params.realref} -I {input.seq} -o {output}"


#rule IndelRealigner:
    #input:
        #bam=resultdir+"{sample}_dedup.bam",
        #ref=hg
        #target=resultdir+"{sample}.interval",
    #output:

    #run:
        #pass

###############################################################################
#                           SINGLE-TIME-RUN RULES                             #
###############################################################################

rule download_reference:
    """download the hg19 human reference genome from 1000genome"""
    output:
        zipped = hg+'.gz',
    version: 0.1
    benchmark:
        "benchmarks/benchmark_downloadreference_subset_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs)
    shell:
        "wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz && "
        "mv human_g1k_v37.fasta.gz {output.zipped} "

rule gunzip_reference:
    input:
        zipped = hg+'.gz'
    output:
        hg
    benchmark:
        "benchmarks/benchmark_gunzip_subset_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs)
    shell:
        "gunzip -k {input.zipped} || true"


rule index_bwa:
    """generate the index of the reference genome for the bwa program"""
    input:
        hg,
    output:
        bwa_indexes,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_index_bwa_subset_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs)
    version: 0.1
    shell:
        "bwa index -a bwtsw {hg}"


rule index_picard:
    """generate the index of the reference genome for the picard program"""
    input:
        hg=hg,
    output:
        hg.replace('fasta', 'dict'),
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_index_picard_subset_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs)
    shell:
        "picard CreateSequenceDictionary R={input.hg} O={output}"


rule index_samtools:
    """generate the index of the reference genome for the samtools and gatk programs"""
    input:
        hg=hg,
    output:
        hg+'.fai'
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_index_samtools_subset_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs)
    shell:
        "samtools faidx {input.hg} "

# https://downloads.sourceforge.net/project/bio-bwa/bwakit/bwakit-0.7.12_x64-linux.tar.bz2
