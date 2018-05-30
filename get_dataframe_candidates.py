import pandas as pd
import itertools
from itertools import *
from Bio import SeqIO
import numpy as np
import argparse
import sys

# This script allows you to merge you two blasy hit files with taxonomic informations with information about distances, GC, coverage etc.

parser = argparse.ArgumentParser(description='Merge and get new dataframe of genes candidates')
parser.add_argument("-c", "--candidates_sp1_sp2", help="introduce de dataframe with informations about dn_ds, coverage etc (output from gff_cov_GC.py)")
parser.add_argument("-b1", "--blast_sp1", help="introduce the blast file with all best hits and taxid inf of the specie1")
parser.add_argument("-b2", "--blast_sp2", help="introduce the blast file with all best hits and taxid ing of the specie2")
parser.add_argument("-s1", "--specie1", help="introduce the name of the specie 1")
parser.add_argument("-s2", "--specie2", help="introduce the name of the specie 2")



args = parser.parse_args()

candidates_sp1_sp2=args.candidates_sp1_sp2
blast_sp1=args.blast_sp1
blast_sp2=args.blast_sp2
specie1=args.specie1 #espèce 1
specie2=args.specie2 #espèce 2

#Ex Usage: 
#python3 get_dataframe_candidates.py -c candidates_genes_pvalue_0035_0042.tab -b1 outfile_0035.tab_top_blast_hits.csv -b2 outfile_0042.tab_top_blast_hits.csv -s1 0035 -s2 0042

augustus_cluster=pd.read_csv(candidates_sp1_sp2,header=0,sep='\t')
blast_sp1=pd.read_table(blast_sp1,header=0,sep=',')
blast_sp2=pd.read_table(blast_sp2,header=0,sep=',')


data = augustus_cluster.merge(blast_sp1, left_on='seq1_id', right_on='qseqid')
data = data.merge(blast_sp2, left_on='seq2_id', right_on='qseqid')


Complet_data=data.drop(columns=["Unnamed: 0", "Mean_length", 'GC_content_seq1','GC_content_seq2','length_x','mismatch_x','gapopen_x','qstart_x','qend_x','sstart_x','send_x',
	'bitscore_x','length_y','mismatch_y','gapopen_y','qstart_y','qend_y','sstart_y','send_y','bitscore_y'])

Complet_data.columns = ['seq1_id','seq2_id','dN','dS','Dist_third_pos','Dist_brute',
 'Length_seq_1','Length_seq_2','scaf_name_seq1','GC_scaff_seq1','cov_depth_seq1','scaf_name_seq2','GC_scaff_seq2','cov_depth_seq2',
 'qseqid_sp1','sseqid_sp1','pident_sp1','evalue_sp1','salltitles_sp1','staxids_sp1','scientific_name_sp1','scomnames_sp1','sskingdoms_sp1','Order_sp1','qseqid_sp2',
 'sseqid_sp2','pident_sp2','evalue_sp2','salltitles_sp2','staxids_sp2','scientific_name_sp2','scomnames_sp2','sskingdoms_sp2','Order_sp2']

#Creat a new dataframe with all informations of these candidates paired sequences:
#Complet_data.to_csv("Complet_candidates_tax.tab",sep='\t')

Summarize_data=Complet_data.drop(columns=['dN','Dist_third_pos','Dist_brute','Length_seq_1','Length_seq_2','scaf_name_seq1','GC_scaff_seq1','cov_depth_seq1','scaf_name_seq2','GC_scaff_seq2','cov_depth_seq2',
	'qseqid_sp1','pident_sp1','staxids_sp1','staxids_sp2','scomnames_sp1','scomnames_sp2','qseqid_sp2','pident_sp2','sskingdoms_sp1','sskingdoms_sp2'])

#Creat a new dataframe with a summary of principals informations of the candidates paired sequences:
#Summarize_data.to_csv("Summarize_candidates_tax.tab",sep='\t')

#This part allows you to count how many candidats genes have been found depending on the trainings parameters:
count_sp1_sp1=0
count_sp1_sp2=0 
count_sp2_sp2=0
count_sp2_sp1=0
count=0


for i in Summarize_data['seq1_id']:
	count+=1
	if '_'+specie1+'_'+specie1 in i:
		count_sp1_sp1 +=1
	if '_'+specie1+'_'+specie2 in i:
		count_sp1_sp2 +=1

for i in Summarize_data['seq2_id']:
	if '_'+specie2+'_'+specie2  in i:
		count_sp2_sp2 +=1
	if '_'+specie2+'_'+specie1  in i:
		count_sp2_sp1 +=1


candidates=pd.read_csv(candidates_sp1_sp2,header=0,sep='\t')
All_candidates = candidates.merge(Complet_data, how='left')
All_candidates.to_csv("All_candidates_with_all_inf.tab",sep='\t')
Summarize_all_candidates = candidates.merge(Summarize_data, how='left')
Summarize_all_candidates.to_csv("Summarize_All_candidates_with_all_inf.tab",sep='\t')


print("Nb_",specie1,"_",specie1," :",count_sp1_sp1,sep="")
print("Nb_",specie1,"_",specie2," :",count_sp1_sp2,sep="")
print("Nb_",specie2,"_",specie2," :",count_sp2_sp2,sep="")
print("Nb_",specie2,"_",specie1," :",count_sp2_sp1,sep="")
print('Nb total :',count_sp1_sp1+count_sp1_sp2+count_sp2_sp2+count_sp2_sp1)

candidates_sp2 = pd.DataFrame(augustus_cluster['seq2_id'])
print(candidates_sp2)


#Allows to only keep sequences after the coverage filter: 
result=All_candidates.merge(candidates_sp2, left_on='seq2_id', right_on='seq2_id')
result.to_csv("All_candidates_cov_all_inf.tab",sep='\t')


