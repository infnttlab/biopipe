import itertools as it
import sys
import os

start = sys.argv[1]
finish = sys.argv[2]

origin_1 = os.path.expanduser("~") + "/data_ref/SRR1611178_1.fastq"
origin_2 = os.path.expanduser("~") + "/data_ref/SRR1611178_2.fastq"
target_1 = "./data/subset_{start}_{finish}_1.fastq".format(start=start,finish=finish)
target_2 = "./data/subset_{start}_{finish}_2.fastq".format(start=start,finish=finish)

files = [(origin_1,target_1),(origin_2,target_2)]

for (origin, target) in files:
    with open(origin, 'r') as origin_file:
        with open(target, 'w') as target_file:
            begin = int(start)*4
            end = int(finish)*4
            partial_lines = it.islice(origin_file, begin, end)
            for line in partial_lines:
                stripped_line = line.strip()
                print(stripped_line, file=target_file)
            
# %%
        
