import itertools as it

start = 0
nreads = 1000000

origin_1 = "SRR1611178_1.fastq"
origin_2 = "SRR1611178_2.fastq"
target_1 = "subset_{start}_{nreads}_1.fastq".format(start=start,nreads=nreads)
target_2 = "subset_{start}_{nreads}_2.fastq".format(start=start,nreads=nreads)

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
        