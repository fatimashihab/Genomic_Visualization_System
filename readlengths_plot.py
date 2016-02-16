import pysam
import matplotlib.pyplot as plt
import os
import sys

def plot(file):
  end=[]
  position=[]
  levels=[-5]
  found=False
  
  
  for read in file.fetch(until_eof=True):
      end.append(read.reference_end)
      position.append(read.pos)
      
      for m,l in enumerate(levels):
          if read.pos > l:
              i = read.reference_end
              levels[m]=i
              found=True
              coeff=1
              break
  
          if not found:
              i= read.reference_end
              levels.append(i)
      
          x=read.pos
          xaxis=[i for i in range (x,i)]
          yaxis=[coeff]*len(xaxis)
          coeff=coeff+2
  
          lines= plt.plot(xaxis,yaxis)
          plt.setp(lines,linewidth=2)
  
  plt.show()

try:
    os.chdir("/home/fatma/Desktop")
    filename=sys.argv[1]
    file= pysam.AlignmentFile(filename, "rb")
    
except IOError:
    try:
        file_path= raw_input ("Please enter file path: \n")
        os.chdir(file_path)

        filename=sys.argv[1]
        file= pysam.AlignmentFile(filename, "rb")

    except OSError:
        print("Path or File Name Error.")
    else:
        position= []
        number_of_reads=[]

        for pileupcolumn in file.pileup(until_eof=True):
            number_of_reads.append(pileupcolumn.n)
            position.append(pileupcolumn.pos)

        #print (number_of_reads)
        #print ("\n %s" %position)

        plt.plot(position, number_of_reads)
        plt.show()

        file.close()


else:
    position= []
    number_of_reads=[]

    for pileupcolumn in file.pileup(until_eof=True):
        number_of_reads.append(pileupcolumn.n)
        position.append(pileupcolumn.pos)

    #print (number_of_reads)
    #print ("\n %s" %position)

    plt.plot(position, number_of_reads)
    plt.show()

    file.close()
