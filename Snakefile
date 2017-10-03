#snakemake --dag | dot -Tpdf > dag.pdf
#snakemake --detailed-summary > provenance.tsv
#snakemake --cores 6 --resources mem=12
shell.executable("/bin/bash")

configfile: "config.yaml"

n_sim = config['n_sim']
cpu_type = config['cpu_type']
thrs = config['threads']
n_cpu = config['n_cpu']

import os
home = os.path.expanduser("~")
hg = home + config['hg']
bwa_indexes = [hg+".bwt", hg+".pac", hg+".amb", hg+".ann", hg+".sa"]
indels_ref = home + config['indels_ref']


datadir = 'data/'
resultdir = 'results/'


samples = [filename for filename in os.listdir('./'+datadir) if filename.endswith('.fastq')]
samples = set("_".join(filename.split('_')[:-1]) for filename in samples)


rule all:
    input:
        expand(resultdir+"{sample}_dedup.bai", sample=samples),
        hg.replace('fasta', 'dict'),
        expand(resultdir+"{sample}_realigned.bam", sample=samples),
    benchmark:
        "benchmarks/benchmark_rule_all_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
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
        "benchmarks/benchmark_mapping_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    threads: 32
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
        "benchmarks/benchmark_sort_picard_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
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
        "benchmarks/benchmark_mark_duplicates_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
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
        "benchmarks/benchmark_build_bam_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
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
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_realigner_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T RealignerTargetCreator -R {params.realref} -I {input.seq} -o {output}"


rule IndelRealigner:
    input:
        bam=resultdir+"{sample}_dedup.bam",
        target=resultdir+"{sample}.interval",
        indels_ref=indels_ref
    output:
        resultdir+"{sample}_realigned.bam",
    params:
        gatk = config['gatk'],    
        #gatk='programs/gatk/GenomeAnalysisTK.jar',
        realref=hg,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_indelrealigner_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T IndelRealigner -R {input.indels_ref} -I {input.bam} -targetIntervals {input.target} -known {input.indels_ref} -o {output}"



###############################################################################
#                           SINGLE-TIME-RUN RULES                             #
###############################################################################

rule download_reference:
    """download the hg19 human reference genome from 1000genome"""
    output:
        zipped = hg+'.gz',
    version: 0.1
    benchmark:
        "benchmarks/benchmark_downloadreference_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz && "
        "mv human_g1k_v37.fasta.gz {output.zipped} "

rule gunzip_reference:
    input:
        zipped = hg+'.gz'
    output:
        hg
    benchmark:
        "benchmarks/benchmark_gunzip_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "gunzip -k {input.zipped} || true"

rule download_indels_ref: 
    """download the indel reference from 1000genome """
    output:
        indel_zipped= indels_ref+'.gz'
    benchmark:
        "benchmarks/benchmark_downloadindels_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        "mv Mills_and_1000G_gold_standard.indels.hg38.vcf.gz {output.indel_zipped}"

rule gunzip_indelref:
    input:
        indel_zipped = indels_ref+'.gz'
    output:
        indels_ref
    benchmark:
        "benchmarks/benchmark_gunzip_indelref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "gunzip -k {input.indel_zipped} || true"

rule index_bwa:
    """generate the index of the reference genome for the bwa program"""
    input:
        hg,
    output:
        bwa_indexes,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_index_bwa_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
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
        "benchmarks/benchmark_index_picard_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
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
        "benchmarks/benchmark_index_samtools_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "samtools faidx {input.hg} "

# https://downloads.sourceforge.net/project/bio-bwa/bwakit/bwakit-0.7.12_x64-linux.tar.bz2
