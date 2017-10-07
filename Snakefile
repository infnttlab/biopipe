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

scripts = '.' + config['scripts']

gatk = home + config['gatk']
#gatk='programs/gatk/GenomeAnalysisTK.jar'

indels_ref = home + config['indels_ref']
dbsnp = home + config['dbsnp']

genotyping_mode = config['hap-caller']['genotyping_mode']
stand_emit_conf = config['hap-caller']['stand_emit_conf']
stand_call_conf = config['hap-caller']['stand_call_conf']

filter_exp_snps = config['hard-filter']['snps']
filter_exp_indels = config['hard-filter']['indels']

convert2annovar = home + config['convert2annovar']
annovar = home + config['annovar']
humandb = home + config['humandb']
build_ver = config['build_ver']
dbsnp_ver = config['dbsnp_ver']
kg_ver = config['kg_ver']
mitochondrial_ver = config['mitochondrial_ver']

datadir = 'data/'
resultdir = 'results/'


samples = [filename for filename in os.listdir('./'+datadir) if filename.endswith('.fastq')]
samples = set("_".join(filename.split('_')[:-1]) for filename in samples)


rule all:
    input:
        #expand(resultdir+"{sample}_recal.bai", sample=samples),
        hg.replace('fasta', 'dict'),
        expand(resultdir+"{sample}.tsv", sample=samples),
    benchmark:
        "benchmarks/benchmark_rule_all_ref_null_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    run:
        pass
        
#########################################################
#                      GATK-LODn                        #
#########################################################

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
        indels_ref=indels_ref
    output:
        resultdir+"{sample}.intervals",
    params:
        gatk = gatk,
        ref=hg,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_realigner_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    threads: 32
    shell:
        "java -jar {params.gatk} -T RealignerTargetCreator -R {params.ref} -I {input.seq} -known {input.indels_ref} -nt {threads} -o {output}"


rule IndelRealigner:
    input:
        bam=resultdir+"{sample}_dedup.bam",
        target=resultdir+"{sample}.intervals",
    output:
        r_bam=resultdir+"{sample}_realigned.bam",
        r_idx=resultdir+"{sample}_realigned.bai",
    params:
        gatk = gatk,
        ref=hg,
        indels_ref=indels_ref
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_indelrealigner_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T IndelRealigner -R {params.ref} -I {input.bam} -targetIntervals {input.target} -known {params.indels_ref} -o {output.r_bam}"

rule BQSR_step_1:
    input:
        r_bam=resultdir+"{sample}_realigned.bam",
        r_idx=resultdir+"{sample}_realigned.bai",
        dbsnp = dbsnp,
    output:
        resultdir+"{sample}_recal_data.table"
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
        "java -jar {params.gatk} -T BaseRecalibrator -R {params.ref} -I {input.r_bam} -knownSites {input.dbsnp} -knownSites {params.indels_ref} -nct {threads} -o {output}"

rule BQSR_step_2:
    input:
        outtable1=resultdir+"{sample}_recal_data.table",
        r_bam=resultdir+"{sample}_realigned.bam",
    output:
        resultdir+"{sample}_post_recal_data.table"
    params:  
        gatk = gatk,
        ref=hg,
        indels_ref=indels_ref,
        dbsnp = dbsnp,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_BQSR2_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    threads: 32
    shell:
        "java -jar {params.gatk} -T BaseRecalibrator -R {params.ref} -I {input.r_bam} -knownSites {params.dbsnp} -knownSites {params.indels_ref} -BQSR {input.outtable1} -nct {threads} -o {output}"

rule BQSR_step_3:
    input:
        outtable1 = resultdir+"{sample}_recal_data.table",
        outtable2 = resultdir+"{sample}_post_recal_data.table",
    output:
        plots = resultdir+"{sample}_recalibrationPlots.pdf",
    params:  
        gatk = gatk,
        ref=hg,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_BQSR3_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "java -jar {params.gatk} -T AnalyzeCovariates -R {params.ref} -before {input.outtable1} -after {input.outtable2} -plots {output.plots}"

rule BQSR_step_4:
    input:
        outtable1 = resultdir+"{sample}_recal_data.table",
        plots = resultdir+"{sample}_recalibrationPlots.pdf",
        r_bam = resultdir+"{sample}_realigned.bam",
    output:
        recal_bam = resultdir+"{sample}_recal.bam",
        recal_bai = resultdir+"{sample}_recal.bai",
    params:  
        gatk = gatk,
        ref=hg,
    conda:
        "envs/config_conda.yaml"
    benchmark:
        "benchmarks/benchmark_BQSR4_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    threads: 32
    shell:
        "java -jar {params.gatk} -T PrintReads -R {params.ref} -I {input.r_bam} -BQSR {input.outtable1} -nct {threads} -o {output.recal_bam}"

#rule BQSR:
#    input:
#        r_bam=resultdir+"{sample}_realigned.bam",
#        r_idx=resultdir+"{sample}_realigned.bai",
#        dbsnp = dbsnp,
#    output:
#        outtable1 = resultdir+"{sample}_recal_data.table",
#        outtable2 = resultdir+"{sample}_post_recal_data.table",
#        plots = resultdir+"{sample}_recalibrationPlots.pdf",
#        recal_bam = resultdir+"{sample}_recal.bam",
#        recal_bai = resultdir+"{sample}_recal.bai",
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
    """GATK Haplotype Caller"""
    input:
        recal_bam = resultdir+"{sample}_recal.bam",
        recal_bai = resultdir+"{sample}_recal.bai",
    output:
        vcf_raw = resultdir+"{sample}_raw.vcf"
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
    """GATK Hard Filtering Variants"""
    ## Extract the snps from the call set
    input:
        vcf_raw = resultdir+"{sample}_raw.vcf"
    output:
        raw_snps = resultdir+"{sample}_snps.vcf"
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
    """GATK Hard Filtering Variants"""
    ## Apply the filter to the SNP call set
    input:
        raw_snps = resultdir+"{sample}_snps.vcf"
    output:
        filt_snps = resultdir+"{sample}_hard_filtered_snps.vcf"
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
    """GATK Hard Filtering Variants"""
    ## Extract the Indels from the call set
    input:
        vcf_raw = resultdir+"{sample}_raw.vcf",
    output:
        raw_indels = resultdir+"{sample}_indels.vcf"
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
    """GATK Hard Filtering Variants"""
    ## Apply the filter to the Indel call set
    input:
        raw_indels = resultdir+"{sample}_indels.vcf",
    output:
        filt_indels = resultdir+"{sample}_hard_filtered_indels.vcf"
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
    """GATK Hard Filtering Variants"""
    ## Combine variants
    input:
        filt_snps = resultdir+"{sample}_hard_filtered_snps.vcf",
        filt_indels = resultdir+"{sample}_hard_filtered_indels.vcf",
    output:
        vcf_filt = resultdir+"{sample}_filtered_variants.vcf"
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
    # convert to annovar format
    input:
        vcf_filt = resultdir+"{sample}_filtered_variants.vcf",
    output:
        outfile = resultdir+"{sample}_filtered_variants.annovar"
    params:  
        convert2annovar = convert2annovar,
        fmt = 'vcf4',
        pars = '',
    benchmark:
        "benchmarks/benchmark_Annotationconvert_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    shell:
        "{params.convert2annovar} -format {params.fmt} {input.vcf_filt} -outfile {output.outfile} -includeinfo '-withzyg' -comment {params.pars}"

###################################################################################
# ./annotate_variation.pl -downdb -buildver hg19 -webfrom annovar snp138 humandb  #
# ./annotate_variation.pl -downdb 1000g2012apr humandb -buildver hg19             #
###################################################################################
    
rule annovar_filter_dbSNP138:
    # annotate
    #dbSNP138
    input:
        outfile = resultdir+"{sample}_filtered_variants.annovar",
    output:
        dbsnp_rmdup = resultdir+"{sample}_rmdup.dbsnp",
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
    # annotate
    #1000g annotation
    input:
        dbsnp_rmdup = resultdir+"{sample}_rmdup.dbsnp",
        outfile = resultdir+"{sample}_filtered_variants.annovar",
    output:
        kg_rmdup = resultdir+"{sample}_rmdup.1000g",
    params:
        annovar = annovar,
        humandb = humandb,
        build_ver = build_ver,
        dbsnp_ver = dbsnp_ver,
        kg_ver = kg_ver,
        mutect = False,
        pars ='-maf 0.05 -reverse',
    benchmark:
        "benchmarks/benchmark_annovarfilter1000g_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    run:
        shell("{params.annovar} -filter {input.outfile}.{params.build_ver}_{params.dbsnp_ver}_filtered -buildver {params.build_ver} -dbtype {params.kg_ver} {params.humandb} {params.pars}")
        suffix = '.%s_%s.sites.\d\d\d\d_\d\d_filtered'%(build_ver,kg_ver[-3:].upper())
        kg_ann = resultdir + [x for x in os.listdir(resultdir) if re.findall(suffix,x)][0]
        shell("awk \'{{print $3,$4,$5,$6,$7,$8,$9,'{params.kg_ver}',$2,$19,$20}}\'"+" {kg_ann}".format(kg_ann=kg_ann)+" > {output.kg_rmdup}")

rule Annotation:
    ## Gene-based annotation
    input:
        dbsnp_rmdup = resultdir+"{sample}_rmdup.dbsnp",
        kg_rmdup = resultdir+"{sample}_rmdup.1000g",
    output: 
        known_file = resultdir+"{sample}_rmdup.known",
        novel_file = resultdir+"{sample}_rmdup.novel",
    params:
        annovar = annovar,
        humandb = humandb,
        build_ver = build_ver,
    benchmark:
        "benchmarks/benchmark_Annotation_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    run:
        suffix = '.%s_%s.sites.\d\d\d\d_\d\d_filtered'%(build_ver,kg_ver[-3:].upper())
        kg_ann = resultdir + [x for x in os.listdir(resultdir) if re.findall(suffix,x)][0]
        shell("cat {input.dbsnp_rmdup} {input.kg_rmdup} > {output.known_file}")
        shell("mv "+"{kg_ann}".format(kg_ann=kg_ann)+" {output.novel_file}")

rule Ann_mitochondrial:
    # annotate
    # Mitochondrial Annotation
    input:
        outfile = resultdir+"{sample}_filtered_variants.annovar",
        known_file = resultdir+"{sample}_rmdup.known",
        novel_file = resultdir+"{sample}_rmdup.novel",
    output:
        k_f = resultdir+"{sample}_rmdup.known.exonic_variant_function", 
        n_f = resultdir+"{sample}_rmdup.novel.exonic_variant_function",
        mit_rmdup = resultdir+"{sample}.mit",
    params:  
        annovar = annovar,
        humandb = humandb,
        build_ver = build_ver,
        mutect = False,
        mitochondrial_ver = mitochondrial_ver,
    benchmark:
        "benchmarks/benchmark_Annmitochondrial_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    run:
        for ann_file in [{input.known_file}, {input.novel_file}]:
            shell("{params.annovar} -geneanno "+"{ann_file}"+" -buildver {params.build_ver} {params.humandb}")
        shell("{params.annovar} -buildver {params.mitochondrial_ver} -dbtype ensGene {input.outfile} {humandb}")
        shell("awk \'{{print $3,$4,$5,$6,$7,$8,$9,'{params.mitochondrial_ver}',$2,$21,$22,$23}}\' {input.outfile}.exonic_variant_function > {output.mit_rmdup}")

rule MakeFinalFile:
    input:
        k_f = resultdir+"{sample}_rmdup.known.exonic_variant_function", 
        n_f = resultdir+"{sample}_rmdup.novel.exonic_variant_function",
        mit_rmdup = resultdir+"{sample}.mit",
    output:
        out = resultdir+"{sample}.tsv",
    params:
        scripts = scripts,
        mutect = False,
        sample_order = ['n','t'],
        dbsnp_freq = True,
        dbsnpFreq = None,
        dbsnpAllele = None,    
    benchmark:
        "benchmarks/benchmark_MakeFinalFile_ref_{sample}" + "_n_sim_{n_sim}_cputype_{cpu_type}_thrs_{thrs}_ncpu_{n_cpu}.txt".format(n_sim=n_sim, cpu_type=cpu_type, thrs=thrs, n_cpu=n_cpu)
    script:
        "{params.scripts}" + "MakeFinalFile.py"
    
        
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
#########################################################
#                   REFERENCE INDEXING                  #
#########################################################

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


#########################################################
#                        THE END                        #
#########################################################