#!/bin/bash
#$ -cwd
#$ -l h_data=64G,h_rt=4:00:00
#$ -N hsc_analysis
#$ -o output.log
#$ -e error.log
#$ -m bea
#$ -M your_email@ucla.edu

# Initialize module system
source /u/local/Modules/default/init/modules.sh

# Load HDF5 and R modules
module load hdf5
module load R/4.2.2-BIO

# Run the analysis
Rscript /u/home/n/neop/hsc_project/mature_expanded_dataset.R

