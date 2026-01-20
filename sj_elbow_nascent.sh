#!/bin/bash
#$ -cwd
#$ -l h_data=128G,h_rt=8:00:00
#$ -N nascent_elbow
#$ -o output_fixed.log
#$ -e error_fixed.log
#$ -m bea
#$ -M neop@g.ucla.edu

# CRITICAL: Initialize the module system first
source /u/local/Modules/default/init/modules.sh

# Now load modules
module load R/4.2.2-BIO
module load hdf5

# Use your personal library
export R_LIBS_USER=~/Rlibs/4.2

# Run the script
Rscript elbow_nascent.R
