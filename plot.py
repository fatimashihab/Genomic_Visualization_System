import pysam
import matplotlib.pyplot as plt

import os
os.chdir("/home/fatma/Desktop/samtools/lib")

file= pysam.AlignmentFile("ex1_sorted.sorted.bam", "rb")

position= []
number_of_reads=[]

for pileupcolumn in file.pileup(until_eof=True):
    number_of_reads.append(pileupcolumn.n)
    position.append(pileupcolumn.pos)

plt.plot(position, number_of_reads)
plt.show()

file.close()
