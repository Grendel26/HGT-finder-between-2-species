#PBS -S /bin/bash
#PBS -l nodes=1:ppn=8:bigmem,mem=30gb
#PBS -e /pandata/bguinet/LEPIWASP/LOGS/Busco_sp1.error
#PBS -o /pandata/bguinet/LEPIWASP/LOGS/Busco_sp1.out
#PBS -N run_busco_sp1
#PBS -q q1day


# Ce script permet de faire tourner le programme Busco (recherche de gènes ultra-conservés 
# dans le génome en une seule copie unique) chez une espèce;
# Ceci permet à la fois d'avoir une idée de la qualité d'assemblage,
# et d'avoir une distribution sous l'hypothèse nulle de divergence entre pairs de séquences
# pour des gènes ayant été transmis par transfert vertical uniquement.

# Fichiers requis: 
# - Génome assemblé sp1 format fasta (assembly)
# - Base de donnée d'espèces d'intérêt (lineage) (à retrouver sur http://busco.ezlab.org/)
# - Dépendances Busco (Augustus, HMM, Blast, Python)

# Prenez soin de définir les noms de variables et les chemins vers les dépendances et fichiers avant de lancer les commandes.

hostname
uname -a

echo debut:
date  #debut 

#Your samp name
SAMP=sp1

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

