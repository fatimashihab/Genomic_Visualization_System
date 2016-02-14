import pysam
import matplotlib.pyplot as plt
import os

def plot(file):
    position= []
    number_of_reads=[]

    for pileupcolumn in file.pileup(until_eof=True):
        number_of_reads.append(pileupcolumn.n)
        position.append(pileupcolumn.pos)

    #print (number_of_reads)
    #print ("\n %s" %position)

    plt.plot(position, number_of_reads)
    plt.ylabel("reads")
    plt.xlabel("contig")
    plt.title("Sam Files Visualization Tool")
    plt.show()

try:
    os.chdir("/home/fatma/Desktop/samtools/lib")
    file= pysam.AlignmentFile("ex1_sorted.sorted.bam", "rb")
except IOError:
    try:
        file_path= raw_input ("Please enter file path: \n")
        os.chdir(file_path)

        file= raw_input("\nPlease enter the bamfile (ex: file_name.bam): \n")
        file=pysam.AlignmentFile(file, "rb")

    except OSError:
        print("Path or File Name Error.")
    else:
        plot(file)
        file.close()


else:
    plot(file)

    file.close()
