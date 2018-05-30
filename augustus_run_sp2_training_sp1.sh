#PBS -S /bin/bash
#PBS -l nodes=1:ppn=8:bigmem,mem=30gb
#PBS -e /pandata/bguinet/LEPIWASP/ACG-sp2_Augustus/LOGS/augustus_training_sp1.error
#PBS -o /pandata/bguinet/LEPIWASP/ACG-sp2_Augustus/LOGS/augustus_training_sp1.out
#PBS -q q1day
#PBS -N augustus_sp2_training_sp1


# Prediction of genes in the genome of species 2 using the sp1 training.

# This script is used to run the Augustus program that will predict genes by an Ab initio approach.
# With Busco run,  training files are available to better predict genes.
# Note that this retraining file (myspecie) must be added to the species file of the Augustus program.

# Usual Busco retraining file path: sp1_busco/run_sp1_BUSCO_v2/augustus_output/retraining_parameters
# The usual path of the file where to transfer the retraining file: / augustus / config / species / myspecie

# Required files:
# - Genome Assembled Fasta Format (assembly)
# - Training file for species 2 (retraining)


# Be sure to define variable names and file paths before issuing commands.


echo debut:
date  #debut 

hostname
uname -a

# 
# 

#Declarations des variables 
SAMP=sp2
PATH=bguinet@pbil-deb

# GENOME (sequence file)
ASSEMBLY=/pandata/varaldi/LEPIWASP/OUT/$SAMP/scaffold.fa

# augustus config path (where the retraining file from busco should be added in the specie file)
PATH=/panhome/bguinet/TOOLS/augustus/config/

#retraining file's name (not a path)
RETRAINING=retraining_sp1
#Augustus programme path
AUGUSTUS=/panhome/bguinet/TOOLS/augustus/bin/augustus


# outputfile's path and name
OUTPUT=/pandata/bguinet/LEPIWASP/ACG-sp2_Augustus/LOGS/run_augustus_sp2_training_sp1.out

$AUGUSTUS --species=$RETRAINING --AUGUSTUS_CONFIG_PATH=$PATH $ASSEMBLY > $OUTPUT


