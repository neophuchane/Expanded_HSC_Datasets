#!/bin/bash
#$ -cwd
#$ -l h_data=64G,h_rt=6:00:00
#$ -N nascent_hsc
#$ -o nascent_output.log
#$ -e nascent_error.log
#$ -m bea
#$ -M neop@g.ucla.edu

# Initialize module system
source /u/local/Modules/default/init/modules.sh

# Load HDF5 and R modules
module load hdf5
module load R/4.2.2-BIO

# Run the analysis
Rscript /u/home/n/neop/hsc_project/expand_nascent.R

