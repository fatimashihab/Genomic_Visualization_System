#!/bin/bash

echo "Batch file to install needed packges for Genome Visualization"

echo "1. Install pip"
sudo apt-get install python-pip
echo ""

echo "2. Install numpy"
sudo pip install numpy
echo ""

echo "3. Install Python-dev"
sudo apt-get install python-dev
echo ""

echo "4. Install biopython"
sudo pip install biopython
echo ""

echo "5. Install docopt"
sudo pip install docopt==0.6.2
echo ""

echo "6. Install pysam"
sudo pip install pysam
echo ""

echo "7. Install plotly"
sudo pip install plotly
sudo pip install plotly --upgrade
echo ""

echo "Right, I'm all done.  Bye bye."
