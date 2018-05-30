#PBS -S /bin/bash
#PBS -l nodes=1:ppn=8:bigmem,mem=30gb
#PBS -e /pandata/bguinet/LEPIWASP/ACG-sp1_Augustus/LOGS/augustus_training_sp2.error
#PBS -o /pandata/bguinet/LEPIWASP/ACG-sp1_Augustus/LOGS/augustus_training_sp2.out
#PBS -q q1day
#PBS -N augustus_sp1_training_sp2


# Prédiction de gènes dans le génome de l'espèce 1 en utilisant l'entrainement de l'espèce 2.

# Ce script permet de faire tourner le programme Augustus qui va prédire des gènes par une approche Ab initio chez une espèce.
# Grâce aux run Busco, les fichiers d'entraînement des espèces sont disponibles
# afin de prédire plus correctement les gènes.
# Notez que ce fichier retraining (myspecie) doit être ajouté dans le fichier species du programme Augustus

# Chemin habituel du fichier retraining de Busco : sp2_busco/run_sp2_BUSCO_v2/augustus_output/retraining_parameters
# Chemin habituel du fichier où transferer le fichier retraining: /augustus/config/species/myspecie

# Fichiers requis: 
# - Génome assemblé format fasta (assembly)
# - Fichier d'entrainement de l'espèce 2 (retraining)


# Prenez soin de définir les noms de variables et les chemins vers les fichiers avant de lancer les commandes.




echo debut:
date  #debut 

hostname
uname -a

# 
# 

#Declarations des variables 
SAMP=sp1
PATH=bguinet@pbil-deb

# GENOME (sequence file)
ASSEMBLY=/pandata/varaldi/LEPIWASP/OUT/$SAMP/scaffold.fa

# augustus config path (where the retraining file from busco should be added in the specie file)
PATH=/panhome/bguinet/TOOLS/augustus/config/

#retraining file's name  (not a path)
RETRAINING=retraining_sp2
#Augustus programme path
AUGUSTUS=/panhome/bguinet/TOOLS/augustus/bin/augustus


# outputfile's path and name
OUTPUT=/pandata/bguinet/LEPIWASP/ACG-sp1_Augustus/LOGS/run_augustus_sp1_training_sp2.out

$AUGUSTUS --species=$RETRAINING --AUGUSTUS_CONFIG_PATH=$PATH $ASSEMBLY > $OUTPUT


