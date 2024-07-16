from Bio import SeqIO
import numpy as np
import sys

arr = []
with open(str(sys.argv[1])) as f:
    for record in SeqIO.parse(f, "fastq"):
        arr.append(len(record.seq))
print(int(np.median(arr)))
