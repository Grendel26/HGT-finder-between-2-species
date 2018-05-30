#PBS -S /bin/bash
#PBS -l nodes=1:ppn=8:bigmem,mem=30gb
#PBS -e /pandata/bguinet/LEPIWASP/LOGS/Busco_sp2.error
#PBS -o /pandata/bguinet/LEPIWASP/LOGS/Busco_sp2.out
#PBS -N run_busco_sp2
#PBS -q q1day


# This script allows you to run the Busco program (search for ultra-conserved genes
# in the genome in one single copy) in one species;
# This allows both an idea of the quality of assembly,
# and to have a distribution under the null hypothesis of divergence between pairs of sequences
# for genes that have been transmitted by vertical transfer only.

# Required files:
# - Assembled genome sp2 format fasta (assembly)
# - Species of Interest Database (lineage) (to be found on http://busco.ezlab.org/)
# - Busco Outbuildings (Augustus, HMM, Blast, Python)

# Be sure to define variable names and paths to dependencies and files before issuing commands.


hostname
uname -a

echo debut:
date  #debut 

#Your samp name
SAMP=sp2

# run BUSCO on genome assembly

#######################
# define path to data #
#######################
# GenomeE (sequence file)
ASSEMBLY=/path/to/your/data/$SAMP/scaffold.fa
# Species data 
LINEAGE=/path/to/your/data/arthropoda_odb9/
# Name output dir
NAME=$SAMP'_BUSCO_v2' # WARNING: do not provide a path.


#########################################
# define PATH to sofwtare used by BUSCO #
#########################################
#Augustus
export PATH=/bin:/usr/bin:/usr/remote/bin:/panhome/varaldi/TOOLS/augustus-3.2.3/
# hmmer
PATH=$PATH:/usr/remote/hmmer-3.1b2/bin
# blast et python
PATH=$PATH:/usr/bin
PATH=$PATH:/panusr/ncbi-blast-2.2.25+/bin/
# augustus
#PATH=$PATH:/panhome/varaldiTOOLS/augustus-3.2.3/bin/ 

export AUGUSTUS_CONFIG_PATH=/panhome/varaldi/TOOLS/augustus-3.2.3/config

#Busco software path
BUSCO=/panhome/varaldi/TOOLS/BUSCO_v2/BUSCO.py

################
# Command line #
################
python $BUSCO -i $ASSEMBLY -o $NAME -l $LINEAGE -m geno -f


echo fin:
date  #fin 

