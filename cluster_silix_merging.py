import pandas as pd
import itertools
from itertools import *
from Bio import SeqIO
import numpy as np
import sys
import argparse



#Script wich allows to deal with Diamond output (.m8 format from silix), Cluster information (.fnodes from silix) and pident mean by Hsp (.net format from silix)
#The next commande lines will only keep cluster with not self-matchs sequences and the best pairs (max pident) within these clusters.

# Print out a message when the program is initiated.
#print('----------------------------------------------------------------\n')
#print('                    Cluster-Silix merging 1.0.\n')
#print('----------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Merge en keep clusters with the best pident')
parser.add_argument("-b", "--blast", help="introduce the blast output (.m8)")
parser.add_argument("-c", "--cluster", help="introduce the cluster file (.fnodes)")
parser.add_argument("-d", "--distance", help="introduce the distance file (.net)")
parser.add_argument("-f1", "--fasta_aa", help="introduce the  fasta file with augustus predicted genes (concatened amino acid)")
parser.add_argument("-f2", "--fasta_dna", help="introduce the  fasta file with augustus predicted genes (concatened dna)")
args = parser.parse_args()


# Variable that stores fasta sequences
blast=args.blast
cluster=args.cluster 
distance=args.distance
fasta_aa=args.fasta_aa
fasta_dna=args.fasta_dna

#Usage xemple:	
#python3 cluster_silix_merging.py -b matches_Augustus_0035_0042.m8 -c cluster_Augustus_0035_0042.fnodes -d matches_Augustus_0035_0042.net -f1 concatenate_0035_0042.faa -f2 concatenate_0035_0042.fna 

blast=pd.read_table(blast,header=None)
blast.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen","qstart", "qend", "sstart", "send", "evalue", "bitscore"]
blast=blast.drop(columns=["mismatch", "gapopen", "evalue", "bitscore"])

#Cluster Dataframe
cluster=pd.read_table(cluster,header=None)
cluster.columns = ["cluster_name", "seq_names"]

#Distance mean dataframe
dist=pd.read_table(distance,header=None)
dist.columns = ["qseqid", "sseqid","pident","coverage"]
dist=dist.drop(columns=["coverage"])

#Including cluster information and distance mean information into one dataframe:
data = cluster.merge(dist, left_on='seq_names', right_on='qseqid')

#Adding for each two remaining dataframe a concatened colomn
data["name_concatened"] = data["qseqid"].map(str) + data["sseqid"]
blast["name_concatened"] = blast["qseqid"].map(str) + blast["sseqid"]
#We do not need these columns anymore
blast=blast.drop(columns=[ "qseqid","sseqid"])

#Including cluster information + distance mean information  + coordinate sequences from blast into one dataframe:
data = data.merge(blast, left_on='name_concatened', right_on='name_concatened')
data=data.drop(columns=["seq_names","name_concatened"])

data.to_csv("dataframe_brute.txt", sep='\t')
#data should only contain pairings of different sequences:
data = data[data['qseqid']!=data["sseqid"]]


#To ignore pairings which have the same substrings in their seqid, the most readable way would be to add data columns with these data:
data['qspec'] = [seqid.split('_')[1] for seqid in data['qseqid'].values]
data['sspec'] = [seqid.split('_')[1] for seqid in data['sseqid'].values]

data_wo_eqSpec = data[data['qspec']!=data['sspec']]


data_wo_eqSpec.to_csv("dataframe.txt", sep='\t')


df=pd.read_table("dataframe.txt",header=0,sep='\t')

#Sorte sequences 
df[['qseqid','sseqid']] = np.sort(df[['qseqid','sseqid']], axis=1)

#droping duplicates
data_wo_eqSpec  = df.drop_duplicates(subset=['qseqid','sseqid'])

#In the end the data should be grouped by cluster-ID and within each group, the maximum of pident is of interest:

data_grpd = data_wo_eqSpec.groupby(['cluster_name'])

data_wo_eqSpec.to_csv("dataframe_not_max.txt",sep='\t')
result=data_wo_eqSpec.loc[data_grpd['pident_x'].idxmax()]


result.to_csv("dataframe_max.txt",sep='\t')


df=pd.read_table("dataframe_max.txt",sep="\t")


output_aa_file = open('clusters1_aa.fasta','w')
output_aa_file2 = open('clusters2_aa.fasta','w')

output_dna_file = open('clusters1_dna.fasta','w')
output_dna_file2 = open('custers2_dna.fasta','w')


#Display into 4 files the dna and amino acid sequences:

#record_dict = SeqIO.to_dict(SeqIO.parse("concatenate_with_busco_names_0035_0042_aa.fa", "fasta"))

qseqid=df.ix[:,3]
sseqid=df.ix[:,4]


record_dict = SeqIO.to_dict(SeqIO.parse(fasta_aa, "fasta"))

for a, b in zip(qseqid,sseqid):
	
	if a in record_dict:
		print(">",record_dict[a].id,file=output_aa_file,sep="")
		print(record_dict[a].seq,file=output_aa_file)
	if a not in record_dict:
		print(a ,"pas dans dic a_aa",b)

	if b in record_dict:
		print(">",record_dict[b].id,file=output_aa_file2,sep="")
		print(record_dict[b].seq,file=output_aa_file2)
	if b not in record_dict:
		print(b ,"pas dans dic b_aa",a)


#Fasta file with all the sequences and their ID
record_dict2 = SeqIO.to_dict(SeqIO.parse(fasta_dna, "fasta"))

qseqid=df.ix[:,3]
sseqid=df.ix[:,4]


for a, b, in zip(qseqid,sseqid):

	if a in record_dict2:
		print(">",record_dict2[a].id,file=output_dna_file)
		print(record_dict2[a].seq,file=output_dna_file)
	if a not in record_dict2:
		print(a ,"pas dans dic a_dna",b)

	if b in record_dict2:
		print(">",record_dict2[b].id,file=output_dna_file2)
		print(record_dict2[b].seq,file=output_dna_file2)
	if b not in record_dict2:
		print(b," pas dans dic b_dna",a)



output_aa_file.close()
output_aa_file.close()
output_aa_file2.close()
output_dna_file.close()
output_dna_file2.close()



