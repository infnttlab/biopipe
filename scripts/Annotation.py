import os
import subprocess
import re

def annovar_filter(annovar,infile, build_ver,dbtype, humandb,outdir='.',
                   err_log=True, pars=''):
    
    if err_log:
        log='2>%s/err.annovar_%s'%(outdir,dbtype)
    else:
        log = ''
    
    subprocess.call(' '.join([annovar, '-filter', infile,'-buildver',
                              build_ver,'-dbtype', dbtype, humandb,
                              pars,log]),shell=True)
    
def process_annovar_out(dbtype, ann, rmdup,mutect=False):
    if mutect:
        command = 'awk \'{print $3,$4,$5,$6,$7,$8,$9,"%s",$2,$17,$18,$19}\''%dbtype
    else:
        command = 'awk \'{print $3,$4,$5,$6,$7,$8,$9,"%s",$2,$19,$20}\''%dbtype
    subprocess.call(' '.join([command, ann,'>', rmdup]),shell=True)

def Annotation(vcf,convert2annovar,annovar,humandb,build_ver, dbsnp_ver,
               kg_ver,mitochondrial_ver,fmt = 'vcf4',outdir='.',name='sample',
               err_log=True,maf=0.05,mutect=False):
    
    if mutect:
        fmt = 'vcf4old'
        pars = '--filter pass'
    else:
        pars = ''
    
    # convert to annovar format
    outdir = os.path.abspath(outdir)
    outfile = '%s/%s.annovar'%(outdir,name)
    
    if err_log:
        log='2>%s/err.convert2annovar'%outdir
    else:
        log = ''
    subprocess.call(' '.join([convert2annovar, '-format',fmt, vcf, '-outfile',
                              outfile, '-includeinfo', '-withzyg',
                              '-comment', pars,log]),shell=True)
    # annotate
    #dbSNP138 and 1000g annotation
    annovar_filter(annovar, outfile,build_ver, dbsnp_ver,humandb,outdir)
    dbsnp_ann = '%s.%s_%s_dropped'%(outfile,build_ver,dbsnp_ver)
    dbsnp_rmdup = '%s/%s_rmdup.dbsnp'%(outdir,name)
    process_annovar_out(dbsnp_ver, dbsnp_ann, dbsnp_rmdup,mutect)
    dbsnp_filt = '%s.%s_%s_filtered'%(outfile,build_ver,dbsnp_ver)
    
    annovar_filter(annovar, dbsnp_filt, build_ver,kg_ver,humandb,outdir,
                   pars='-maf 0.05 -reverse')
    suffix = '.%s_%s.sites.\d\d\d\d_\d\d_filtered'%(build_ver,
                                                    kg_ver[-3:].upper())
    
    kg_ann = [x for x in os.listdir(outdir) if re.findall(suffix,x)][0]
    kg_ann = outdir+'/'+kg_ann
    kg_rmdup = '%s/%s_rmdup.1000g'%(outdir,name)
    process_annovar_out(kg_ver, kg_ann, kg_rmdup,mutect)
    
    known_file = '%s/%s_rmdup.known'%(outdir,name)
    novel_file = '%s/%s_rmdup.novel'%(outdir,name)
    
    os.system('cat %s %s > %s'%(dbsnp_rmdup,kg_rmdup,known_file))
    
    #Turning your life more easy =) -> one changing the name of novel variants file
    os.system('mv %s %s'%(kg_ann,novel_file))
    
    ## Gene-based annotation
    for ann_file in [known_file, novel_file]:
        subprocess.call(' '.join([annovar,'-geneanno',ann_file,'-buildver',
                              build_ver,humandb]), shell=True)
    
    #Mitochondrial Annotation
    
    subprocess.call(' '.join([annovar, '-buildver', mitochondrial_ver,
                              '-dbtype ensGene', outfile, humandb,
                              log]), shell=True)
    if mutect:
        command = 'awk \'{print $3,$4,$5,$6,$7,$8,$9,"%s",$2,$19,$20,$21}\''
        command = command%mitochondrial_ver
    else:
        command = 'awk \'{print $3,$4,$5,$6,$7,$8,$9,"%s",$2,$21,$22,$23}\''
        command = command%mitochondrial_ver
    mit_ann = '%s.exonic_variant_function'%outfile
    mit_rmdup = '%s/%s.mit'%(outdir,name)
    subprocess.call(' '.join([command, mit_ann, '>',mit_rmdup]),shell=True)
    
    ann_files = ['%s.exonic_variant_function'%known_file,
                 '%s.exonic_variant_function'%novel_file,
                 mit_rmdup]
    
    return(ann_files)

Annotation(snakemake.input['vcf'],snakemake.params['convert2annovar'],snakemake.params['annovar'],snakemake.params['humandb'],snakemake.params['build_ver'],
           snakemake.params['dbsnp_ver'], snakemake.params['kg_ver'],snakemake.params['mitochondrial_ver'],snakemake.params['fmt'],snakemake.params['outdir'],
           snakemake.params['name'], snakemake.params['err_log'],snakemake.params['maf'],snakemake.params['mutect'])