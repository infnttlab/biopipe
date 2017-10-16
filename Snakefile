#snakemake --dag | dot -Tpdf > dag.pdf
#snakemake --detailed-summary > provenance.tsv
#snakemake --cores 6 --resources mem=12

shell.executable("/bin/bash") # The pipeline works only on /bin/bash

configfile: "config.yaml" # Set configuration file

# Extract the main directory
import os
home = os.path.expanduser("~")

#####################################
#  Set parameters from config file  #
#####################################

# Label's parameters
n_sim = config['label_params']['n_sim']
cpu_type = config['label_params']['cpu_type']
thrs = config['label_params']['threads']
n_cpu = config['label_params']['n_cpu']
 
# References
hg = home + config['ref-files']['hg'] # Human Genome Reference
bwa_indexes = [hg+".bwt", hg+".pac", hg+".amb", hg+".ann", hg+".sa"] # Indexes for hg
indels_ref = home + config['ref-files']['indels_ref'] # Set of known indels
dbsnp = home + config['ref-files']['dbsnp'] # SNP database
cosmic = home + config['ref-files']['cosmic'] # Catalog of somatic mutation in cancer
humandb = home + config['ref-files']['humandb'] # Annovar human databases folder
build_ver = config['ref-files']['build_ver'] # Set build version
dbsnp_ver = config['ref-files']['dbsnp_ver'] # Set SNP database version
kg_ver = config['ref-files']['kg_ver'] # Set  version # Set ethnicity groups based on time
mitochondrial_ver = config['ref-files']['mitochondrial_ver'] # Set parameter command for annotating mitochondria variants

# Folders
scripts = config['folders']['scripts']
datadir = config['folders']['datadir']
normal_dir = config['folders']['datadir_normal']
tumour_dir = config['folders']['datadir_tumour']
resultdir = config['folders']['resultdir']

# Softwares
gatk = home + config['softwares']['gatk']
muTect = home + config['softwares']['muTect']
annovar = home + config['softwares']['annovar']
convert2annovar = home + config['softwares']['convert2annovar']

# Sample details
platform = config['sample-details']['platform'] # Set platform for mapping
library = config['sample-details']['library'] # Set library for mapping
target = home + config['sample-details']['target'] # Set target intervals for exome analysis

# HaplotypeCaller parameters
genotyping_mode = config['hap-caller']['genotyping_mode'] # Set genotyping mode
stand_emit_conf = config['hap-caller']['stand_emit_conf'] # stand_emit_conf 10.0 means that it won’t report any potential SNPs with a quality below 10.0; 
stand_call_conf = config['hap-caller']['stand_call_conf'] # but unless they meet the quality threshold set by -stand_call_conf (30.0, in this case), they will be listed as failing the quality filter

# Hard Filter parameters
filter_exp_snps = config['hard-filter']['snps'] # Set command to define thresholds for snps
filter_exp_indels = config['hard-filter']['indels'] # Set command to define thresholds for indels

# Annovar databases
annovar_dbs = [home + config['annovar_dbs']['hg19_db'] , home + config['annovar_dbs']['snp138_db']]

#LODn parameters
min_n_cov = config['LODn']['min_n_cov'] #Set the minimum read depth for normal and tumor samples
min_t_cov = config['LODn']['min_t_cov']

#######################
#   Get fastq files   #
#######################

patients = [folder for folder in os.listdir(datadir)]
p_folders = [normal_dir, tumour_dir]
SAMPLES = [filename for Dir in os.listdir(datadir) for sub_dir in p_folders if os.path.isdir(datadir+Dir+sub_dir) for filename in os.listdir(datadir+Dir+sub_dir) if filename.endswith('.fastq')]
SAMPLES = set("_".join(filename.split('_')[:-1]) for filename in SAMPLES)
sicks = [folder for folder in os.listdir(datadir) if (os.path.isdir(datadir+folder+normal_dir)) and (os.path.isdir(datadir+folder+tumour_dir))] 

wildcard_constraints:
    sample = "("+"|".join(SAMPLES)+")",
    sick = "("+"|".join(sicks)+")",

samples_directory = {}
sample_to_patient = {}
for patient in patients:
    for sub_dir in p_folders:
        if os.path.isdir(datadir+patient+sub_dir):
            for filename in os.listdir(datadir+patient+sub_dir):
               if filename.endswith('.fastq'):
                   filename = "_".join(filename.split('_')[:-1])
                   samples_directory[filename] = datadir+patient+sub_dir
                   sample_to_patient[filename] = patient 

def set_dir1(wildcards):
    return (samples_directory[wildcards]+wildcards+'_1.fastq')
    
def set_dir2(wildcards):
    return (samples_directory[wildcards]+wildcards+'_2.fastq')


sicks_list = {}
for sick in sicks:
    sicks_list[sick] = {}
    for sub_dir in p_folders:
        name = [filename for filename in os.listdir(datadir+sick+sub_dir) if filename.endswith('.fastq')][0]       
        name = "_".join(name.split('_')[:-1])
        sicks_list[sick][sub_dir] = name
   

def get_bam(wildcards, sample_type='N'):    
    if sample_type=='N': 
        return (resultdir+sicks_list[wildcards][normal_dir]+"_recal.bam")
    elif sample_type=='T':   
        return (resultdir+sicks_list[wildcards][tumour_dir]+"_recal.bam")    

        
def get_code(wildcards):
    return (sicks_list[wildcards][normal_dir])
        
def get_lodn_infile(wildcards,format):
    if format == 'table':
        return (resultdir+sicks_list[wildcards][tumour_dir]+".tsv")
    elif format == 'vcf':
        return (resultdir+sicks_list[wildcards][tumour_dir]+"_filtered_variants.vcf")

#######################
# Workflow final rule #
#######################

rule all:
    input:
        expand(resultdir+"{sample}"+".tsv", sample=SAMPLES),
        expand(resultdir+ 'mutect_ann/' + "{sick}"+".tsv", sick=sicks),
        expand(resultdir+"{sick}"+"_lodn_vcf.tsv",sick=sicks),
        expand(resultdir+"{sick}"+"_lodn_table.tsv",sick=sicks),
    benchmark:
        "benchmarks/benchmark_rule_all_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    run:
        pass
        
#########################################################
#                      GATK-LODn                        #
#########################################################

#############################
#        gatk_pipe          #
#############################


rule mapping:
    """
    This tool maps the samples to the reference genome.
    It does append a header necessary for the GATK analysis.
    """
    input:
        bwa_indexes,
        #sample1 = "{sample}"+"_1.fastq",
        #sample2 = "{sample}"+"_2.fastq",
        sample1 = lambda wildcards:set_dir1(wildcards.sample),
        sample2 = lambda wildcards:set_dir2(wildcards.sample),
    output:
        outfile = resultdir+"{sample}"+".sam",
    params:
        reference = hg,
        library = library,
        platform = platform,
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
    """
    This tool sorts the input SAM or BAM file by coordinate.
    """
    input:
        r = resultdir+"{sample}"+".sam",
    output:
        resultdir+"{sample}"+"_sorted.bam",
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_sort_picard_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard SortSam INPUT={input.r} OUTPUT={output}.tmp SORT_ORDER=coordinate"
        " && [ -s {output}.tmp ] && mv {output}.tmp {output}"


rule mark_duplicates:
    """
    This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. 
    """
    input:
        r=resultdir+"{sample}"+"_sorted.bam",
    output:
        resultdir+"{sample}"+"_dedup.bam",
    params:
        outdir=resultdir,
        metricsfile=resultdir+"{sample}"+"_dedup.metrics.txt",
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
    """
    This tool creates an index file for the input BAM that allows fast look-up of data in a BAM file, like an index on a database.
    """
    input:
        r=resultdir+"{sample}"+"_dedup.bam",
    output:
        resultdir+"{sample}"+"_dedup.bai",
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_build_bam_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "picard BuildBamIndex INPUT={input.r} OUTPUT={output}"


rule realigner_target_creator:
    """
    This tool defines intervals to target for local realignment.
    """
    input:
        hg.replace('fasta', 'dict'),
        idx=resultdir+"{sample}"+"_dedup.bai",
        ref=hg+'.fai',
        indels_ref=indels_ref,
        gatk = gatk,
    output:
        resultdir+"{sample}"+".intervals",
    params:
        seq=resultdir+"{sample}"+"_dedup.bam",
        ref=hg,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_realigner_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    threads: 32
    shell:
        "java -jar {input.gatk} -T RealignerTargetCreator -R {params.ref} -I {params.seq} -known {input.indels_ref} -nt {threads} -o {output}"


rule IndelRealigner:
    """
    This tool performs local realignment of reads around indels.
    """
    input:
        target=resultdir+"{sample}"+".intervals",
    output:
        r_bam=resultdir+"{sample}"+"_realigned.bam",
        r_idx=resultdir+"{sample}"+"_realigned.bai",
    params:
        gatk = gatk,
        ref=hg,
        indels_ref=indels_ref,
        bam=resultdir+"{sample}"+"_dedup.bam",
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_indelrealigner_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T IndelRealigner -R {params.ref} -I {params.bam} -targetIntervals {input.target} -known {params.indels_ref} -o {output.r_bam}"

rule BQSR_step_1:
    """
    This tool detects systematic errors in base quality scores.
    This step produces a recalibrated data table.
    """
    input:
        r_bam=resultdir+"{sample}"+"_realigned.bam",
        r_idx=resultdir+"{sample}"+"_realigned.bai",
        dbsnp = dbsnp,
    output:
        outtable1 = resultdir+"{sample}"+"_recal_data.table"
    params:  
        gatk = gatk,
        ref=hg,
        indels_ref=indels_ref
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_BQSR1_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    threads: 32
    shell:
        "java -jar {params.gatk} -T BaseRecalibrator -R {params.ref} -I {input.r_bam} -knownSites {input.dbsnp} -knownSites {params.indels_ref} -nct {threads} -o {output.outtable1}"

rule BQSR_step_2:
    """
    This tool detects systematic errors in base quality scores.
    This step produces a post recalibrated data table.
    """
    input:
        outtable1 = resultdir+"{sample}"+"_recal_data.table",
    output:
        outtable2 = resultdir+"{sample}"+"_post_recal_data.table"
    params:  
        gatk = gatk,
        ref=hg,
        indels_ref=indels_ref,
        dbsnp = dbsnp,
        r_bam=resultdir+"{sample}"+"_realigned.bam",
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_BQSR2_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    threads: 32
    shell:
        "java -jar {params.gatk} -T BaseRecalibrator -R {params.ref} -I {params.r_bam} -knownSites {params.dbsnp} -knownSites {params.indels_ref} -BQSR {input.outtable1} -nct {threads} -o {output.outtable2}"

rule BQSR_step_3:
    """
    This tool creates plots to visualize base recalibration results. 
    """
    input:
        outtable2 = resultdir+"{sample}"+"_post_recal_data.table",
    output:
        plots = resultdir+"{sample}"+"_recalibrationPlots.pdf",
    params:  
        gatk = gatk,
        ref=hg,
        outtable1 = resultdir+"{sample}"+"_recal_data.table",
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_BQSR3_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T AnalyzeCovariates -R {params.ref} -before {params.outtable1} -after {input.outtable2} -plots {output.plots}"

rule BQSR_step_4:
    """
    This tool writes out sequence read data.
    """
    input:
        plots = resultdir+"{sample}"+"_recalibrationPlots.pdf",
    output:
        recal_bam = resultdir+"{sample}"+"_recal.bam",
        recal_bai = resultdir+"{sample}"+"_recal.bai",
    params:  
        gatk = gatk,
        ref=hg,
        r_bam = resultdir+"{sample}"+"_realigned.bam",
        outtable1 = resultdir+"{sample}"+"_recal_data.table",
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_BQSR4_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    threads: 32
    shell:
        "java -jar {params.gatk} -T PrintReads -R {params.ref} -I {params.r_bam} -BQSR {params.outtable1} -nct {threads} -o {output.recal_bam}"

#rule BQSR:
#    input:
#        r_bam=resultdir+"{sample}"+"_realigned.bam",
#        r_idx=resultdir+"{sample}"+"_realigned.bai",
#        dbsnp = dbsnp,
#    output:
#        outtable1 = resultdir+"{sample}"+".recal_data.table",
#        outtable2 = resultdir+"{sample}"+".post_recal_data.table",
#        plots = resultdir+"{sample}"+"_recalibrationPlots.pdf",
#        recal_bam = resultdir+"{sample}"+"_recal.bam",
#        recal_bai = resultdir+"{sample}"+"_recal.bai",
#    params:  
#        gatk = gatk,
#        ref=hg,
#        indels_ref=indels_ref
#    conda:
#        "envs/config_conda.yaml"
#    benchmark:
#        "benchmarks/benchmark_BQSR_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
#    threads: 32
#    shell:
#        "java -jar {params.gatk} -T BaseRecalibrator -R {params.ref} -I {input.r_bam} -knownSites {input.dbsnp} -knownSites {params.indels_ref} -nct {threads} -o {output.outtable1} && "
#        "java -jar {params.gatk} -T BaseRecalibrator -R {params.ref} -I {input.r_bam} -knownSites {input.dbsnp} -knownSites {params.indels_ref} -BQSR {output.outtable1} -nct {threads} -o {output.outtable2} && "
#        "java -jar {params.gatk} -T AnalyzeCovariates -R {params.ref} -before {output.outtable1} -after {output.outtable2} -plots {output.plots} && "
#        "java -jar {params.gatk} -T PrintReads -R {params.ref} -I {input.r_bam} -BQSR {output.outtable1} -nct {threads} -o {output.recal_bam}"

rule HaplotypeCaller:
    """
    GATK Hard Filtering Variants
    This tool calls germline SNPs and indels via local re-assembly of haplotypes.
    """
    input:
        recal_bam = resultdir+"{sample}"+"_recal.bam",
        recal_bai = resultdir+"{sample}"+"_recal.bai",
    output:
        vcf_raw = resultdir+"{sample}"+"_raw.vcf"
    params:  
        gatk = gatk,
        ref=hg,
        genotyping_mode = genotyping_mode,
        stand_emit_conf = stand_emit_conf,
        stand_call_conf = stand_call_conf,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_HaplotypeCaller_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T HaplotypeCaller -R {params.ref} -I {input.recal_bam} --genotyping_mode {params.genotyping_mode} -stand_emit_conf {params.stand_emit_conf} -stand_call_conf {params.stand_call_conf} -o {output.vcf_raw}" 
        
# Since GATK 3.7 -stand_emit_conf {params.stand_emit_conf} has been deprecated

rule HardFilter_1SNPs:
    """
    GATK Hard Filtering Variants
    This tool extracts the snps from the call set.
    """
    input:
        vcf_raw = resultdir+"{sample}"+"_raw.vcf"
    output:
        raw_snps = resultdir+"{sample}"+"_snps.vcf"
    params:  
        gatk = gatk,
        ref=hg,
        selectType = '-selectType SNP'
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_HardFilter1SNPs_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T SelectVariants -R {params.ref} -V {input.vcf_raw} {params.selectType} -o {output.raw_snps}"     

rule HardFilter_2SNPs:
    """
    GATK Hard Filtering Variants
    This tool applies the filter to the SNP call set.
    """
    input:
        raw_snps = resultdir+"{sample}"+"_snps.vcf"
    output:
        filt_snps = resultdir+"{sample}"+"_hard_filtered.snps.vcf"
    params:  
        gatk = gatk,
        ref=hg,
        filter_exp_snps = filter_exp_snps,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_HardFilter2SNPs_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T VariantFiltration -R {params.ref} -V {input.raw_snps} {params.filter_exp_snps} -o {output.filt_snps}"

rule HardFilter_1Indels:
    """
    GATK Hard Filtering Variants
    This tool extracts the Indels from the call set.
    """
    input:
        vcf_raw = resultdir+"{sample}"+"_raw.vcf",
    output:
        raw_indels = resultdir+"{sample}"+"_indels.vcf"
    params:  
        gatk = gatk,
        ref=hg,
        selectType = '-selectType INDEL',
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_HardFilter1Indels_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T SelectVariants -R {params.ref} -V {input.vcf_raw} {params.selectType} -o {output.raw_indels}"

rule HardFilter_2Indels:
    """
    GATK Hard Filtering Variants
    This tool applies the filter to the Indel call set.
    """
    input:
        raw_indels = resultdir+"{sample}"+"_indels.vcf",
    output:
        filt_indels = resultdir+"{sample}"+"_hard_filtered.indels.vcf"
    params:  
        gatk = gatk,
        ref=hg,
        filter_exp_indels = filter_exp_indels,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_HardFilter2Indels_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T VariantFiltration -R {params.ref} -V {input.raw_indels} {params.filter_exp_indels} -o {output.filt_indels}"

rule HardFilter_Combine:
    """
    GATK Hard Filtering Variants
    This tool combines variants.
    """
    input:
        filt_snps = resultdir+"{sample}"+"_hard_filtered.snps.vcf",
        filt_indels = resultdir+"{sample}"+"_hard_filtered.indels.vcf",
    output:
        vcf_filt = resultdir+"{sample}"+"_filtered_variants.vcf"
    params:  
        gatk = gatk,
        ref=hg,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_HardFilterCombine_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T CombineVariants -R {params.ref} --variant:snps {input.filt_snps} --variant:indels {input.filt_indels} -o {output.vcf_filt} -genotypeMergeOptions PRIORITIZE -priority snps,indels"

rule Annotation_convert:
    """
    This step converts to annovar format the vcf files.
    """
    input:
        vcf_filt = resultdir+"{sample}"+"_filtered_variants.vcf",
        annovar_dbs = annovar_dbs,
    output:
        outfile = resultdir+"{sample}"+"_filtered_variants.annovar"
    params:  
        convert2annovar = convert2annovar,
        fmt = 'vcf4',
        pars = '',
    benchmark:
        "benchmarks/benchmark_Annotationconvert_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "{params.convert2annovar} -format {params.fmt} {input.vcf_filt} -outfile {output.outfile} -includeinfo '-withzyg' -comment {params.pars}"

    
rule annovar_filter_dbSNP138:
    """
    This step annotate variants based on dbSNP138.
    """
    input:
        outfile = resultdir+"{sample}"+"_filtered_variants.annovar",
    output:
        dbsnp_rmdup = resultdir+"{sample}"+"_rmdup.dbsnp",
    params:  
        annovar = annovar,
        humandb = humandb,
        build_ver = build_ver,
        dbsnp_ver = dbsnp_ver,
        mutect = False,
        pars = ''
    benchmark:
        "benchmarks/benchmark_annovarfilterdbSNP138_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "{params.annovar} -filter {input.outfile} -buildver {params.build_ver} -dbtype {params.dbsnp_ver} {params.humandb} {params.pars} && "
        "awk \'{{print $3,$4,$5,$6,$7,$8,$9,'{params.dbsnp_ver}',$2,$19,$20}}\' {input.outfile}.{params.build_ver}_{params.dbsnp_ver}_dropped > {output.dbsnp_rmdup}"


rule annovar_filter_1000g:
    """
    This step annotate variants based on 1000g database.
    """
    input:
        dbsnp_rmdup = resultdir+"{sample}"+"_rmdup.dbsnp",
    output:
        kg_rmdup = resultdir+"{sample}"+"_rmdup.1000g",
    params:
        annovar = annovar,
        humandb = humandb,
        build_ver = build_ver,
        dbsnp_ver = dbsnp_ver,
        kg_ver = kg_ver,
        mutect = False,
        pars ='-maf 0.05 -reverse',
        outfile = resultdir+"{sample}"+"_filtered_variants.annovar",
    benchmark:
        "benchmarks/benchmark_annovarfilter1000g_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    run:
        shell("{params.annovar} -filter {params.outfile}.{params.build_ver}_{params.dbsnp_ver}_filtered -buildver {params.build_ver} -dbtype {params.kg_ver} {params.humandb} {params.pars}")
        suffix = '.%s_%s.sites.\d\d\d\d_\d\d_filtered'%(build_ver,kg_ver[-3:].upper())
        kg_ann = resultdir + [x for x in os.listdir(resultdir) if re.findall(suffix,x)][0]
        shell("awk \'{{print $3,$4,$5,$6,$7,$8,$9,'{params.kg_ver}',$2,$19,$20}}\'"+" {kg_ann}".format(kg_ann=kg_ann)+" > {output.kg_rmdup}")

rule Gene_annotation:
    """
    Gene-based annotation.
    """
    input:
        kg_rmdup = resultdir+"{sample}"+"_rmdup.1000g",
    output: 
        known_file = resultdir+"{sample}"+"_rmdup.known",
        novel_file = resultdir+"{sample}"+"_rmdup.novel",
    params:
        annovar = annovar,
        humandb = humandb,
        build_ver = build_ver,
        dbsnp_rmdup = resultdir+"{sample}"+"_rmdup.dbsnp",
    benchmark:
        "benchmarks/benchmark_Geneann_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    run:
        suffix = '.%s_%s.sites.\d\d\d\d_\d\d_filtered'%(build_ver,kg_ver[-3:].upper())
        kg_ann = resultdir + [x for x in os.listdir(resultdir) if re.findall(suffix,x)][0]
        shell("cat {params.dbsnp_rmdup} {input.kg_rmdup} > {output.known_file}")
        shell("mv "+"{kg_ann}".format(kg_ann=kg_ann)+" {output.novel_file}")

rule Ann_mitochondrial:
    """
    Mitochondrial Annotation.
    """
    input:
        known_file = resultdir+"{sample}"+"_rmdup.known",
        novel_file = resultdir+"{sample}"+"_rmdup.novel",
    output:
        k_f = resultdir+"{sample}"+"_rmdup.known.exonic_variant_function", 
        n_f = resultdir+"{sample}"+"_rmdup.novel.exonic_variant_function",
        mit_rmdup = resultdir+"{sample}"+".mit",
    params:  
        annovar = annovar,
        humandb = humandb,
        build_ver = build_ver,
        mutect = False,
        mitochondrial_ver = mitochondrial_ver,
        outfile = resultdir+"{sample}"+"_filtered_variants.annovar",
    benchmark:
        "benchmarks/benchmark_Annmitochondrial_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    run:
        for ann_file in [{input.known_file}, {input.novel_file}]:
            shell("{params.annovar} -geneanno "+"{ann_file}"+" -buildver {params.build_ver} {params.humandb}")
        shell("{params.annovar} -buildver {params.mitochondrial_ver} -dbtype ensGene {params.outfile} {humandb}")
        shell("awk \'{{print $3,$4,$5,$6,$7,$8,$9,'{params.mitochondrial_ver}',$2,$21,$22,$23}}\' {params.outfile}.exonic_variant_function > {output.mit_rmdup}")

rule MakeFinalFile:
    """
    This step makes the final file.
    """
    input:
        k_f = resultdir+"{sample}"+"_rmdup.known.exonic_variant_function", 
        n_f = resultdir+"{sample}"+"_rmdup.novel.exonic_variant_function",
        mit_rmdup = resultdir+"{sample}"+".mit",
    output:
        out = resultdir+"{sample}"+".tsv",
    params:
        scripts = scripts,
        mutect = False,
        sample_order = ['n','t'],
        dbsnp_freq = True,
        dbsnpFreq = None,
        dbsnpAllele = None,
        code = None, # Not necessary except for muTect
        vcf = None,  # Not necessary except for muTect  
    benchmark:
        "benchmarks/benchmark_MakeFinalFile_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    script:
        "{params.scripts}" + "MakeFinalFile.py"


##############################
#        mutect_pipe         #
##############################



rule muTect:
    input:
        muTect = muTect,
        cosmic = cosmic,
        fixed_target = target + "_fixed.bed",
        normal_bam = lambda wildcards: get_bam(wildcards.sick,'N'),
        tumour_bam = lambda wildcards: get_bam(wildcards.sick,'T'),
    output:
        vcf = resultdir+"{sick}"+"_mutect.vcf",
        coverage_out = resultdir+"{sick}"+"_coverage.wig",
    params:
        dbsnp = dbsnp,
        ref = hg,
    conda:
        "envs/config_conda_muTect.yaml"
    benchmark:
        "benchmarks/benchmark_muTect_ref_{sick}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -Xmx50g -jar {input.muTect} --analysis_type MuTect --reference_sequence {params.ref} --cosmic {input.cosmic} --intervals {input.fixed_target} --input_file:normal {input.normal_bam} --input_file:tumor {input.tumour_bam} --vcf {output.vcf} --coverage_file {output.coverage_out}"

rule Annotation_convert_mutect:
    """
    This step converts to annovar format the vcf files.
    """
    input:
        vcf = resultdir+"{sick}"+"_mutect.vcf",
        annovar_dbs = annovar_dbs,
    output:
        outfile = resultdir+ 'mutect_ann/' + "{sick}"+".annovar"
    params:  
        convert2annovar = convert2annovar,
        fmt = 'vcf4old',
        pars = '--filter pass',
    benchmark:
        "benchmarks/benchmark_AnnotationconvertM_ref_{sick}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "{params.convert2annovar} -format {params.fmt} {input.vcf} -outfile {output.outfile} -includeinfo '-withzyg' -comment {params.pars}"
    
rule annovar_filter_dbSNP138_mutect:
    """
    This step annotate variants based on dbSNP138.
    """
    input:
        outfile = resultdir+ 'mutect_ann/' + "{sick}"+".annovar",
    output:
        dbsnp_rmdup = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.dbsnp",
    params:  
        annovar = annovar,
        humandb = humandb,
        build_ver = build_ver,
        dbsnp_ver = dbsnp_ver,
        mutect = True,
        pars = ''
    benchmark:
        "benchmarks/benchmark_annovarfilterdbSNP138M_ref_{sick}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "{params.annovar} -filter {input.outfile} -buildver {params.build_ver} -dbtype {params.dbsnp_ver} {params.humandb} {params.pars} && "
        "awk \'{{print $3,$4,$5,$6,$7,$8,$9,'{params.dbsnp_ver}',$2,$17,$18,$19}}\' {input.outfile}.{params.build_ver}_{params.dbsnp_ver}_dropped > {output.dbsnp_rmdup}"

rule annovar_filter_1000g_mutect:
    """
    This step annotate variants based on 1000g database.
    """
    input:
        dbsnp_rmdup = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.dbsnp",
    output:
        kg_rmdup = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.1000g",
    params:
        annovar = annovar,
        humandb = humandb,
        build_ver = build_ver,
        dbsnp_ver = dbsnp_ver,
        kg_ver = kg_ver,
        mutect = True,
        pars ='-maf 0.05 -reverse',
        outfile = resultdir+ 'mutect_ann/' + "{sick}"+".annovar",
    benchmark:
        "benchmarks/benchmark_annovarfilter1000gM_ref_{sick}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    run:
        shell("{params.annovar} -filter {params.outfile}.{params.build_ver}_{params.dbsnp_ver}_filtered -buildver {params.build_ver} -dbtype {params.kg_ver} {params.humandb} {params.pars}")
        suffix = '.%s_%s.sites.\d\d\d\d_\d\d_filtered'%(build_ver,kg_ver[-3:].upper())
        kg_ann = resultdir + 'mutect_ann/' + [x for x in os.listdir(resultdir + 'mutect_ann') if re.findall(suffix,x)][0]
        shell("awk \'{{print $3,$4,$5,$6,$7,$8,$9,'{params.kg_ver}',$2,$17,$18,$19}}\'"+" {kg_ann}".format(kg_ann=kg_ann)+" > {output.kg_rmdup}")

rule Gene_annotation_mutect:
    """
    Gene-based annotation.
    """
    input:
        kg_rmdup = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.1000g",
    output: 
        known_file = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.known",
        novel_file = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.novel",
    params:
        annovar = annovar,
        humandb = humandb,
        build_ver = build_ver,
        dbsnp_rmdup = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.dbsnp",
    benchmark:
        "benchmarks/benchmark_GeneannM_ref_{sick}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    run:
        suffix = '.%s_%s.sites.\d\d\d\d_\d\d_filtered'%(build_ver,kg_ver[-3:].upper())
        kg_ann = resultdir + 'mutect_ann/' + [x for x in os.listdir(resultdir + 'mutect_ann/') if re.findall(suffix,x)][0]
        shell("cat {params.dbsnp_rmdup} {input.kg_rmdup} > {output.known_file}")
        shell("mv "+"{kg_ann}".format(kg_ann=kg_ann)+" {output.novel_file}")

rule Ann_mitochondrial_mutect:
    """
    Mitochondrial Annotation.
    """
    input:
        known_file = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.known",
        novel_file = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.novel",
    output:
        k_f = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.known.exonic_variant_function", 
        n_f = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.novel.exonic_variant_function",
        mit_rmdup = resultdir+ 'mutect_ann/' + "{sick}"+".mit",
    params:  
        annovar = annovar,
        humandb = humandb,
        build_ver = build_ver,
        mutect = True,
        mitochondrial_ver = mitochondrial_ver,
        outfile = resultdir+ 'mutect_ann/' + "{sick}"+".annovar",
    benchmark:
        "benchmarks/benchmark_AnnmitochondrialM_ref_{sick}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    run:
        for ann_file in [{input.known_file}, {input.novel_file}]:
            shell("{params.annovar} -geneanno "+"{ann_file}"+" -buildver {params.build_ver} {params.humandb}")
        shell("{params.annovar} -buildver {params.mitochondrial_ver} -dbtype ensGene {params.outfile} {humandb}")
        shell("awk \'{{print $3,$4,$5,$6,$7,$8,$9,'{params.mitochondrial_ver}',$2,$19,$20,$21}}\' {params.outfile}.exonic_variant_function > {output.mit_rmdup}")

rule MakeFinalFile_mutect:
    """
    This step makes the final file.
    """
    input:
        k_f = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.known.exonic_variant_function", 
        n_f = resultdir+ 'mutect_ann/' + "{sick}"+"_rmdup.novel.exonic_variant_function",
        mit_rmdup = resultdir+ 'mutect_ann/' + "{sick}"+".mit",
    output:
        out = resultdir+ 'mutect_ann/' + "{sick}"+".tsv",
    params:
        scripts = scripts,
        mutect = True,
        sample_order = ['n','t'], # sample_order can change for muTect
        dbsnp_freq = True,
        dbsnpFreq = None,
        dbsnpAllele = None, 
        code = lambda wildcards: get_code(wildcards.sick),
        vcf = resultdir+"{sick}"+"_mutect.vcf",   
    benchmark:
        "benchmarks/benchmark_MakeFinalFileM_ref_{sick}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    script:
        "{params.scripts}" + "MakeFinalFile.py"
   
###########################
#        run_lodn         #
###########################

rule lodn_table:
    """
    This step applies the LODn classifier on table.
    """
    input:
        infile = lambda wildcards: get_lodn_infile(wildcards.sick,'table'),
        normal_bam = lambda wildcards: get_bam(wildcards.sick,'N'),
    output:
        outfile = resultdir+"{sick}"+"_lodn_table.tsv",
    params:
        fmt = 'table',
        scripts = scripts,
        min_n_cov = min_n_cov,
        min_t_cov = min_t_cov,
    conda:
        "envs/config_conda.yaml"   
    benchmark:
        "benchmarks/benchmark_LODntable_ref_{sick}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    script:
        "{params.scripts}" + "LODn_table.py"

rule lodn_vcf:
    """
    This step applies the LODn classifier on vcf.
    """
    input:
        infile = lambda wildcards: get_lodn_infile(wildcards.sick,'vcf'),
        normal_bam = lambda wildcards: get_bam(wildcards.sick,'N'),
    output:
        outfile = resultdir+"{sick}"+"_lodn_vcf.tsv",
    params:
        fmt = 'vcf',
        scripts = scripts,
        min_n_cov = min_n_cov,
        min_t_cov = min_t_cov,   
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_LODnvcf_ref_{sick}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    script:
        "{params.scripts}" + "LODn_vcf.py"

#rule LODn:
#    """
#    This step applies the LODn classifier.
#    """
#    input:
#        infile = ,
#        normal = ,
#    output:
#        outfile = ,
#    params:
#        fmt = config['fmt'],
#        min_n_cov = ,
#        min_t_cov = ,   
#    benchmark:
#        "benchmarks/benchmark_LODn_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
#    script:
#        "{params.scripts}" + "LODn.py"
        
###############################################################################
#                           SINGLE-TIME-RUN RULES                             #
###############################################################################

#########################################################
#                   DOWNLOADING FILES                   #
#########################################################

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
        "wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz && "
        "mv Mills_and_1000G_gold_standard.indels.b37.vcf.gz {output.indel_zipped}"

rule gunzip_indelref:
    input:
        indel_zipped = indels_ref+'.gz'
    output:
        indels_ref
    benchmark:
        "benchmarks/benchmark_gunzip_indelref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "gunzip -k {input.indel_zipped} || true"
        
        
rule download_dbsnp: 
    """download dbsnp_138.b37 from 1000genome """
    output:
        dbsnp_zipped= dbsnp+'.gz'
    benchmark:
        "benchmarks/benchmark_downloaddbsnp_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz && "
        "mv dbsnp_138.b37.vcf.gz {output.dbsnp_zipped}"

rule gunzip_dbsnp:
    input:
        dbsnp_zipped= dbsnp+'.gz'
    output:
        dbsnp
    benchmark:
        "benchmarks/benchmark_gunzip_dbsnp_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "gunzip -k {input.dbsnp_zipped} || true"


# https://downloads.sourceforge.net/project/bio-bwa/bwakit/bwakit-0.7.12_x64-linux.tar.bz2

rule download_cosmic:
    """download  cosmic from broadinstitute"""
    output:
        cosmic=cosmic,
    shell:
        "wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/b37_cosmic_v54_120711.vcf && "
        "mv b37_cosmic_v54_120711.vcf {output.cosmic}"

rule download_target:
    """download target from illumina"""
    output:
        target = target + ".bed",
    shell:
        "wget https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_expandedexome_targetedregions.bed && "
        "mv nexterarapidcapture_expandedexome_targetedregions.bed {output.target}"

rule fix_target:
    """fix target: remove prefix chr in first column and last column(name)"""
    input:
        target = target + ".bed",
    output:
        fixed_target = target + "_fixed.bed",
    shell:
        "sed 's/^chr//' {input.target} > no_chr.bed &&"
        "awk \'{{print $1,$2,$3}}\' no_chr.bed > {output.fixed_target} &&"
        "rm no_chr.bed"
        

rule download_annovar_databases:
    input:
        annovar = annovar,
    output:
        annovar_dbs,
    params:
        humandb = humandb,
    shell:
        "{input.annovar} -downdb -buildver hg19 -webfrom annovar snp138 {params.humandb} && "
        "{input.annovar} -downdb 1000g2012apr {params.humandb} -buildver hg19"
  
        
#########################################################
#                   REFERENCE INDEXING                  #
#########################################################

rule index_bwa:
    """
    Generate the index of the reference genome for the bwa program.
    """
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
    """
    Generate the index of the reference genome for the picard program.
    """
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
    """
    Generate the index of the reference genome for the samtools and gatk programs.
    """
    input:
        hg=hg,
    output:
        hg+'.fai',
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_index_samtools_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "samtools faidx {input.hg} "

#########################################################
#                   CHECKING REQUIREMENTS               #
#########################################################

rule check_GATK:
    output:
        gatk,
    priority: 4
    shell:
        "echo 'Error. Genome Analysis ToolKit not found in softwares directory.' && "
        "exit 1"

rule check_muTect:
    output:
        muTect,
    priority: 3
    shell:
        "echo 'Error. muTect not found in softwares directory.' && "
        "exit 1"

rule check_Annovar:
    output:
        annovar,
    priority: 2
    shell:
        "echo 'Error. Annovar not found in softwares directory.' && "
        "exit 1"
        
#########################################################
#                        THE END                        #
#########################################################