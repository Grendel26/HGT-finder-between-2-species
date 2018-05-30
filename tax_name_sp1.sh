#PBS -S /bin/bash
#PBS -l nodes=1:ppn=8
#PBS -e /pandata/bguinet/LEPIWASP/blast_database/tax_name_candidates_sp1.error
#PBS -o /pandata/bguinet/LEPIWASP/blast_database/tax_name_candidates_sp1.out
#PBS -q q1day
#PBS -N tax_name_candidates_sp1

#This script allows to get taxid informations in the output blast from diamond for the specie #

source /panhome/bguinet/miniconda3/bin/activate #activate your env
export PYTHONPATH=$PYTHONPATH:/panhome/bguinet/miniconda3/lib/python3.6/site-packages #Get all packages needed
diamond_tab_output=/pandata/bguinet/LEPIWASP/blast_database/matches_sp1_candidates.m8  #Your output blast of your candidates genes 
Diamond_blast_to_taxid=/pandata/bguinet/LEPIWASP/blast_database/public_scripts-master/Diamond_BLAST_add_taxonomic_info/Diamond_blast_to_taxid2.py #path to the python programm 

taxid=/pandata/bguinet/LEPIWASP/blast_database/prot.accession2taxid 

categories=/pandata/bguinet/LEPIWASP/blast_database/categories.dmp

names=/pandata/bguinet/LEPIWASP/blast_database/names.dmp

description=/pandata/bguinet/LEPIWASP/blast_database/acc_to_des.tab

$Diamond_blast_to_taxid -i $diamond_tab_output -t $taxid -c $categories -n $names -d $description -o /pandata/bguinet/LEPIWASP/blast_database/outfile_sp1.tab

