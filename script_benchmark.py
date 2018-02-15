import pandas as pd  
import os
import re
import sys

#benchmark_{name}_subset_{sample}_n_sim_{n_sim}_cputype_{cpu_type}_Totthrs_{thrs}_Rulethrs_{Rthrs}_ncpu{n_cpu}.txt

tablename = sys.argv[1]

list_of_df = []
for file in os.listdir():
    if file==tablename:
        df = pd.read_csv(file, sep='\t')
        list_of_df.append(df)
    
list_of_benchmark = []
for file in os.listdir('./'+'benchmarks'):
    if file.startswith("benchmark_") and file.endswith(".txt"):
        list_of_benchmark.append(file) 
        
for benchmark in list_of_benchmark:
    df = pd.read_csv('./benchmarks/'+benchmark, sep='\t')
    regex = "benchmark_(\w+)_ref_(\w+)_n_sim_(\w+)_cputype_(\w+)_Totthrs_(\w+)_Rulethrs_(\w+)_ncpu_(\w+).txt"
    benchmark_fixed = benchmark.replace("-","_")
    info = re.findall(regex, benchmark_fixed)[0]
    print(info)    
    df['rule'] = info[0]
    if info[1]!="null":
        df['Reference subject'] = info[1]
    df['n_sim'] = info[2]
    df['cpu_type'] = info[3]
    df['Snakemake_threads'] = info[4]
    df['Rule_threads'] = info[5]
    df['#cpu'] = info[6]
    list_of_df.append(df)
x = pd.concat(list_of_df)
x.to_csv(tablename,sep='\t',index = False)
y = pd.read_csv(tablename,sep='\t')
y = y.sort_values(['n_sim'],ascending=True)
y.to_csv(tablename,sep='\t',index = False)
# %%
