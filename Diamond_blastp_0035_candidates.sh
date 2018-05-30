#PBS -S /bin/bash
#PBS -l nodes=1:ppn=8
#PBS -e /pandata/bguinet/LEPIWASP/blast_database/Diamond_blastp_sp1_candidate.error
#PBS -o /pandata/bguinet/LEPIWASP/blast_database/Diamond_blastp_sp1_candidate.out
#PBS -q q1day
#PBS -N diamond_blastp


diamond=/panhome/bguinet/miniconda3/bin/diamond 
nr=/panbanques/diamond/nr/nr.dmnd
candidates_aa_sp1=/pandata/bguinet/LEPIWASP/blast_database/candidates_aa_pvalue_0035.fasta

$diamond blastp -d $nr -q $candidates_aa_sp1 -o /pandata/bguinet/LEPIWASP/blast_database/matches_0035_candidates.m8

