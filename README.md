#Visualization-System-for-Large-Genomic-Datasets

Genome Visualization is Python-based tool that visualizes or plots reads in SAM/BAM files to aid in the analysis of data through different
plotting manners with features enhancing analyzing effects & to aid in the efficacy of reading and manipulating files in the SAM/BAM format.

Genome Visualization

System Requirements:
--------------------------------
Genome Visualization requires Python (2.7 or greater) and Cython (0.22 or greater).
It has not been tested on many other platforms.

Genome Visualization Tool is supported on the following operating systems:
Linux 32-bit (x86) and 64-bit (x86_x64)
   
Prerequisities:
--------------------------------
Genome Visualization is executed primarily depending on Python packages dealing with sam files and plotting techiniques.
Main Packages used are: Pysam, Biopython, Plotly, Docopt, and Numpy.

To install the needed python package, run the batch file setup.sh in your terminal by setting the path to file's location, and running it in your terminal as follows:

$ ./setup.sh 

Usage Instructions:
--------------------------------
Usage:

vst.py (-view FILE) [(-ref NAME --start VALUE --end VALUE )]  [-reff FASTAFILE] [--mp]
vst.py (-h | --help)

Visualizes or plots reads in SAM/BAM files to aid in the analysis of data through different tracks with features enhancing analyzing effects.


Arguments:
  FILE       input file to be visualized [BAM format]
  NAME       Name of the reference chromosome
  VALUE      start and end regions' values
  FASTAFILE  input fasta reference genome


Options:
  -view FILE                   Imports and Views BAM files
  -ref  NAME                   Name of the reference chromosome
  -s,--start VALUE             Determines specific starting region in file
  -e, --end VALUE              Determines the ending region in file (optional)[default: end]
  -reff FASTAFILE              Plots reference genome
  --mp                          Displays matepairs among plots
  -h, --help                   Shows help document and quit
 


*	*	 *	 *	 *	 *	 *	*	*	*	*	



In the terminal use the Usage commands with the 'python' command beforehand.

The '-view' option is obligatory followed by the file intended to be visualized with its full path. The '-view' option opens "Bam" file formats.

You can run the command with just the file, and it will be fully plotted, or you may choose to plot or visualize the file from&to certain specified regions. Visualization of reads' length is displayed in the resulted image. Black reads are correct seq. complement reads, while green reads are reverse complemented. 

To access the specific regions in a file, use the '-ref', '--start' and '--end' options. Once you use one of them, the other will be obligtory.
Use the '-s' or '--start' option followed by the position of the read you would like to start with.
Use the '-e' or '--end' option followed by the position of the read you would like to end at.
And the '-ref' option followed by the name of the reference genome that reads align to.

To plot the reference genome track, use the '-reff' option followed by the file path to the reference genome in "FASTA" format.

To use the MatePairs feature, make sure to type '--mp' in your commandline before you run it. No arguments required.

Use the '-h' or '--help' option after vst.py in your terminal to view the proper usage of VST.



## Author: Fatima Hassan Shihab 
## Feedback: fatima.shihab29@gmail.com
