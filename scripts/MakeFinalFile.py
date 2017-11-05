import os
import subprocess
import pandas as pd
from pandas import DataFrame
from collections import defaultdict

def get_order_vcf (normal_name, infile):
    
    for line in open(infile,'r').readlines():
        if line.startswith('#CHROM'):
            order = line.rstrip().split('\t')[-2:]
            if order.index(normal_name) == 1:
                return(['t','n'])
            if order.index(normal_name) == 0:
                return(['n','t'])

def parse_known (fi_known, mutect=False,sample_order=['n','t']):
    genotype = {'0/0':'hom','0/1': 'het', '1/1': 'hom'}
    
    cmd = "cat %s | tr [:blank:] '\t'"%fi_known
    tmp = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE,shell=True)
    (output, err) = tmp.communicate()
    output = output.decode('utf8')
    known = output.split('\n')
    if known[-1] == '':
        known = known[:-1]
    dic = defaultdict(dict)
    if mutect:
        fmt = 'GT:AD:BQ:DP:FA:SS'
        cols = ['line', 'type','category','annotation','chr', 'pos.start',
                'pos.end','ref.base','alt.base','genotype','QUAL','database',
                'name_var', 'format','info_%s'%sample_order[0],
                'info_%s'%sample_order[1]]
    else:
        fmt = 'GT:AD:DP:GQ:PL'
        cols = ['line', 'type','category','annotation','chr', 'pos.start',
                'pos.end','ref.base','alt.base','genotype','QUAL','database',
                'name_var', 'format','info']
    for i in range(len(known)-1):
        vec = known[i].split('\t')
        for k in range(len(vec)):
            dic[i][cols[k]] = vec[k]
    dt = DataFrame.from_dict(dic, orient='index')
    dt = dt[dt['format'] == fmt]
    if mutect:
        rf = {'t':[],'n':[]}
        at = {'t':[],'n':[]}
        gt = [x.split(':')[0] for x in dt['info_t']]
        ## ugly solution
        ## bug: different genotype coding found in vcf. The code was set for
        ## genotypes coded as in the genotype dictionary
        try:
            dt['genotype'] = [genotype[x] for x in gt]
        except:
            dt['genotype'] = float('nan')

        for i in ['t', 'n']:
            rf[i] = [x.split(':')[1].split(',')[0] for x in dt['info_%s'%i]]
            at[i] = [x.split(':')[1].split(',')[1] for x in dt['info_%s'%i]]
            dt['%s_cov.ref'%i] = rf[i]
            dt['%s_cov.alt'%i] = at[i]               
    else:
        dt['cov.ref'] = [x.split(':')[1].split(',')[0] for x in dt['info']]
        dt['cov.alt'] = [x.split(':')[1].split(',')[1] for x in dt['info']]
    dt['known.flag'] = 1
    return(dt)

def parse_novel (infile,mutect=False,sample_order=['n','t']):
    genotype = {'0/0':'hom','0/1': 'het', '1/1': 'hom'}
    cmd = "cat %s | tr [:blank:] '\t'"%infile
    tmp = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE, shell=True)
    (output, err) = tmp.communicate()
    output = output.decode('utf8')
    data = output.split('\n')
    if data[-1] == '':
        data = data[:-1]
    dic = defaultdict(dict)
    if mutect:
        fmt = 'GT:AD:BQ:DP:FA:SS'
        cols = ['line', 'type','category','annotation','chr', 'pos.start',
                'pos.end','ref.base','alt.base','genotype','QUAL','database',
                'name_var', 'format','info_%s'%sample_order[0],
                'info_%s'%sample_order[1]]
        for i in range(len(data)):
            vec = data[i].split('\t')
            for k in range(len(vec[:11])):
                dic[i][cols[k]] = vec[k]
            dic[i]['format'] = vec[-3]
            dic[i]['info_n'] = vec[-2]
            dic[i]['info_t'] = vec[-1]
    else:
        fmt = 'GT:AD:DP:GQ:PL'
        cols = ['line', 'type','category','annotation','chr', 'pos.start',
                'pos.end','ref.base','alt.base','genotype','QUAL','database',
                'name_var', 'format','info']
        for i in range(len(data)-1):
            vec = data[i].split('\t')
            for k in range(len(vec[:11])):
                dic[i][cols[k]] = vec[k]
            dic[i]['format'] = vec[-2]
            dic[i]['info'] = vec[-1]
    dt = DataFrame.from_dict(dic, orient='index')
    dt = dt[dt['format'] == fmt]
    if mutect:
        rf = {'t':[],'n':[]}
        at = {'t':[],'n':[]}
        gt = [x.split(':')[0] for x in dt['info_t']]
        ## ugly solution
        ## bug: different genotype coding found in vcf. The code was set for
        ## genotypes coded as in the genotype dictionary
        try:
            dt['genotype'] = [genotype[x] for x in gt]
        except:
            dt['genotype'] = float('nan')
        for i in ['t', 'n']:
            rf[i] = [x.split(':')[1].split(',')[0] for x in dt['info_%s'%i]]
            at[i] = [x.split(':')[1].split(',')[1] for x in dt['info_%s'%i]]
            dt['%s_cov.ref'%i] = rf[i]
            dt['%s_cov.alt'%i] = at[i]        
    else:
        dt['cov.ref'] = [x.split(':')[1].split(',')[0] for x in dt['info']]
        dt['cov.alt'] = [x.split(':')[1].split(',')[1] for x in dt['info']]
    dt['name_var'] = 'novel'
    dt['known.flag'] = 0
    return(dt)
    
   
def parse_mit (fi_known,mutect=False,sample_order=['n','t']):  
    genotype = {'0/0':'hom','0/1': 'het', '1/1': 'hom'}
    cmd = "cat %s | tr [:blank:] '\t'"%fi_known
    tmp = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE,shell=True)
    (output, err) = tmp.communicate()
    output = output.decode('utf8')
    data = output.split('\n')
    ## check if last line is empty
    if data[-1] == '':
        data = data[:-1]
    dic = defaultdict(dict)
    if mutect:
        fmt = 'GT:AD:BQ:DP:FA:SS'
        cols = ['category','annotation','chr', 'pos.start',
                'pos.end','ref.base','alt.base','database','type',
                'format','info_%s'%sample_order[0], 
                'info_%s'%sample_order[1]]
        for i in range(len(data)):
            vec = data[i].rstrip().split('\t')
            #ixs = [2,3,4,5,6,7,8,10,1,18,19,20]
            ixs = [0,1,2,3,4,5,6,7,8,9,10,11]
            for k in range(len(ixs)):
                dic[i][cols[k]] = vec[ixs[k]]
    else:
        fmt = 'GT:AD:DP:GQ:PL'
        cols = ['category','annotation','chr', 'pos.start',
                'pos.end','ref.base','alt.base','database','type',
                'format','info']
        for i in range(len(data)):
            vec = data[i].rstrip().split('\t')
            for k in range(len(vec)):
                dic[i][cols[k]] = vec[k]
    dt = DataFrame.from_dict(dic, orient='index')
    dt = dt[dt['format'] == fmt]
    if mutect:
        rf = {'t':[],'n':[]}
        at = {'t':[],'n':[]}
        for i in ['t', 'n']:
            rf[i] = [x.split(':')[1].split(',')[0] for x in dt['info_%s'%i]]
            at[i] = [x.split(':')[1].split(',')[1] for x in dt['info_%s'%i]]
            dt['%s_cov.ref'%i] = rf[i]
            dt['%s_cov.alt'%i] = at[i]        
    else:
        dt['cov.ref'] = [x.split(':')[1].split(',')[0] for x in dt['info']]
        dt['cov.alt'] = [x.split(':')[1].split(',')[1] for x in dt['info']]
    dt['genotype'] = 'Mit'
    dt['name_var'] = 'Mit'
    dt['known.flag'] = 0
    return(dt)

def MakeFinalFile (fi_known,fi_novel,fi_mit,outfile='sample.tsv',mutect=False,
                   sample_order=['n','t'],dbsnp_freq=True, dbsnpFreq = None,
                   dbsnpAllele=None, code = None, vcf = None):
    
    ## defining column names
    cols = []
    if mutect:
        cols = ['name_var','type','chr', 'pos.start', 'pos.end','ref.base',
                'alt.base','genotype','annotation','t_cov.ref','t_cov.alt',
                'n_cov.ref','n_cov.alt','known.flag']
    else:
        cols = ['name_var','type','chr', 'pos.start', 'pos.end','ref.base',
                'alt.base','genotype','annotation','cov.ref','cov.alt',
                'known.flag']
    if mutect:
        sample_order = get_order_vcf (code, vcf)
        
    test_size = 0
    for i in [fi_novel, fi_known, fi_mit]:
        test_size += os.path.getsize(i)
    if test_size == 0:
        final = DataFrame(columns=cols)
    else:
        dt = {}
        if os.path.getsize(fi_novel) > 0:
            dt['novel'] = parse_novel(fi_novel,mutect=mutect,
                                      sample_order=sample_order)
            dt['novel']['known.flag'] = 0
        if os.path.getsize(fi_known) > 0:
            dt['known'] = parse_known(fi_known,mutect=mutect,
                                      sample_order=sample_order)
            dt['known']['known.flag'] = 1
        if os.path.getsize(fi_mit) > 0:
            dt['mit'] = parse_mit(fi_mit,mutect=mutect,
                                  sample_order=sample_order)
            dt['mit']['known.flag'] = 0
        for i in dt.keys():
            if mutect:
                dt[i]['genotype'] = 'unknown'
            dt[i] = dt[i][cols]
        
        final = pd.concat(dt.values(), ignore_index=True)
        
        if dbsnp_freq:
            if dbsnpFreq != None and dbsnpAllele != None:
                ## prepare dict
                dbsnp = defaultdict(dict)
                for i in open(dbsnpFreq).readlines():
                    vec = i.rstrip().split('\t')
                    dbsnp[vec[0]]['allele_id'] = vec[1]
                    dbsnp[vec[0]]['freq'] = vec[3]
                    dbsnp[vec[0]]['allele'] = ''
                    dbsnp[vec[0]]['allele_rev'] = ''
                alleles = defaultdict(dict)
                for i in open(dbsnpAllele).readlines():
                    vec = i.rstrip().split('\t')
                    alleles[vec[0]]['allele'] = vec[1]
                    alleles[vec[0]]['allele_rev'] = vec[3]
                for rs in dbsnp.keys():
                    try:
                        id1 = dbsnp[rs]['allele_id']
                        dbsnp[rs]['allele'] = alleles[id1]['allele']
                        id2 = alleles[dbsnp[rs]['allele_id']]['allele_rev']
                        dbsnp[rs]['allele_rev'] = alleles[id2]['allele']
                    except:
                        next
                        
            ## check freq
            final['dbsnp.freq'] = float('NaN')
            for i in final.index:
                var = final['name_var'].loc[i]
                if var.startswith('rs'):
                    try:
                        final['dbsnp.freq'].loc[i] = dbsnp_dic[var[2:]]
                    except:
                        next
        
    final.to_csv(outfile, sep='\t',header=True, index=None)



MakeFinalFile (snakemake.input['k_f'],snakemake.input['n_f'],snakemake.input['mit_rmdup'],snakemake.output['out'],snakemake.params['mutect'],
                   snakemake.params['sample_order'],snakemake.params['dbsnp_freq'], snakemake.params['dbsnpFreq'],
                   snakemake.params['dbsnpAllele'], snakemake.params['code'],snakemake.params['vcf'])