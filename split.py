import itertools as it

n_times = 5

start = 0
range_width = 100000

for i in range(1, n_times+1):
    nreads = i * range_width
    finish = nreads + start
    origin_1 = "~/data_ref/SRR1611178_1.fastq"
    origin_2 = "~/data_ref/SRR1611178_2.fastq"
    target_1 = "~/biopipe/data/subset_{start}_{finish}_1.fastq".format(start=start,finish=finish)
    target_2 = "~/biopipe/data/subset_{start}_{finish}_2.fastq".format(start=start,finish=finish)
    
    files = [(origin_1,target_1),(origin_2,target_2)]
    
    for origin, target in files:
        with open(origin, 'r') as origin_file:
            with open(target, 'w') as target_file:
                size = 4
                begin = start*size
                end = (start+nreads)*size
                partial_lines = it.islice(origin_file, begin, end)
                for line in partial_lines:
                    stripped_line = line.strip()
                    print(stripped_line, file=target_file)
            
# %%
        