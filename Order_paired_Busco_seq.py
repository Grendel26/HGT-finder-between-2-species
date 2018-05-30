import pandas as pd
import itertools
from itertools import *
from Bio import SeqIO
# Module to parse arguments
import argparse
import os, shutil
from shutil import copyfile

import fileinput
import glob


#Allows from the two Busco output files, to recover 4 files containing only the complete Busco genes in common between the 2 species

# -contains the concatenated aa and dna sequences of sp1 'concatened_'+specie1'
# -contains the concatenated aa and dna sequences of sp1  concatened_'+specie2'
# -fasta dna file sp1 "specie1+".faa"
# -fasta dna file sp2 "pecie1+".fna"
# -fasta aa file sp1 "specie2+".faa"
# -fasta aa file sp2 "specie2+".fna"


# Print out a message when the program is initiated.
#print('----------------------------------------------------------------------\n')
#print('                    Classifier of paired sequences from BUSCO 1.0.\n')
#print('----------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow to classifie paired sequence in order')
parser.add_argument("-s1", "--specie1", help="introduce the desired name of the specie 1, ex : 0035")
parser.add_argument("-s2", "--specie2", help="introduce the desired name of the specie 2, ex : 0042")
parser.add_argument("-d1", "--directory1", help="introduce the directory where are present the Busco files of sp1")
parser.add_argument("-d2", "--directory2", help="introduce the directory where are present the Busco files of sp2")
parser.add_argument("-t1", "--table1", help="introduce the BUSCO summary table  of sp1")
parser.add_argument("-t2", "--table2", help="introduce the BUSCO summary table  of sp2")
#Exemple usage: python3 Order_paired_seq.py  -s1 0035 -s2 0042 -d1 /Users/etudiant/Desktop/Horizon_project/Buscou.out/single_copy_busco_sequences_0035/ -d2 /Users/etudiant/Desktop/Horizon_project/Buscou.out/single_copy_busco_sequences_0042/

args = parser.parse_args()
# Variable that stores fasta sequences
specie1=args.specie1 
specie2=args.specie2 
directory1=args.directory1
directory2=args.directory2
table1=args.table1
table2=args.table2

#First, creat a new folders with all complete Busco sequences:
dir_path = os.path.dirname(os.path.realpath(__file__))
new_directory_sp1= dir_path+"/new_directory_"+specie1+"/"
new_directory_sp2= dir_path +"/new_directory_"+specie2+"/"
if not os.path.exists(new_directory_sp1):
    os.makedirs(new_directory_sp1)
if not os.path.exists(new_directory_sp2):
    os.makedirs(new_directory_sp2)


#Only keep into a dataframe comptes Busco genes 
df1=pd.read_table(table1,sep="\t")
df1=df1[df1["Status"].str.contains("Fragmented")==False]
df1=df1[df1['Status'].str.contains("Duplicated")==False]
df1=df1[df1['Status'].str.contains("Missing")==False]
Number_complet_sp1=df1['Status'].count()

df2=pd.read_table(table2,sep="\t")
df2=df2[df2["Status"].str.contains("Fragmented")==False]
df2=df2[df2['Status'].str.contains("Duplicated")==False]
df2=df2[df2['Status'].str.contains("Missing")==False]
Number_complet_sp2=df2['Status'].count()


#Convert the Busco_id into a list
liste1=[]
for i in df1["Busco id"]:
	liste1.append(i)

liste2=[]
for i in df2["Busco id"]:
	liste2.append(i)

#count=0
#scan the directory 
for filename in os.listdir(directory1):
#match file radix against a set (extracted from dataframe or whatever) for fast lookup
    if os.path.splitext(filename)[0] in liste1:
       #count+=1
       # move into a new directory 
       shutil.copy(directory1 + filename,new_directory_sp1)

#count=0
#scan the directory 
for filename in os.listdir(directory2):
#match file radix against a set (extracted from dataframe or whatever) for fast lookup
    if os.path.splitext(filename)[0] in liste2:
       #count+=1
       # move into a new directory 
       shutil.copy(directory2 + filename,new_directory_sp2)


os.system('cat '+new_directory_sp1+'* > concatened_'+specie1)
os.system('cat '+new_directory_sp2+'* > concatened_'+specie2)

#------------------------------------
Dna_aa_1_species_concatened="concatened_"+specie1
Dna_aa_2_species_concatened="concatened_"+specie2

output_aa_sp1 = open(specie1+".faa",'w')
output_dna_sp1 = open(specie1+".fna",'w')
output_aa_sp2 = open(specie2+".faa",'w')
output_dna_sp2 = open(specie2+".fna",'w')

#Display into 4 files the dna and amino acid sequences:

liste1=[]
for record in SeqIO.parse(Dna_aa_1_species_concatened, "fasta"):
	liste1.append(record)

liste2=liste1[0::2]#only aa sequences
liste3=liste1[1::2]#only dna sequences

SeqIO.write(liste2,specie1+".faa","fasta")#Creat a fasta file with only aa seqs
SeqIO.write(liste3,specie1+".fna","fasta")#Creat a fasta file with only dna seqs

liste3=[]
for record in SeqIO.parse(Dna_aa_2_species_concatened, "fasta"):
	liste3.append(record)
liste5=liste3[0::2]#only aa sequences
liste6=liste3[1::2]#only dna sequences

SeqIO.write(liste5,specie2+".faa","fasta")#Creat a fasta file with only aa seqs
SeqIO.write(liste6,specie2+".fna","fasta")#Creat a fasta file with only dna seqs

sp1_aa = SeqIO.parse(specie1+".faa", "fasta")
sp1_dna = SeqIO.parse(specie1+".fna", "fasta")
p2_aa = SeqIO.parse(specie2+".faa", "fasta")
sp2_dna = SeqIO.parse(specie2+".fna", "fasta")


#Ne garder que le nom de séquence sans le lien de télachargement 
with open(specie1+".faa") as original, open(specie1+".faa2", 'w') as corrected:
    records = SeqIO.parse(specie1+".faa", 'fasta')
    for record in records:
    	record.description = record.id = record.id.split(':', 1)[0]
    	SeqIO.write(record, corrected, 'fasta')

with open(specie1+".fna") as original, open(specie1+".fna2", 'w') as corrected:
    records = SeqIO.parse(specie1+".fna", 'fasta')
    for record in records:
    	record.description = record.id = record.id.split(':', 1)[0]
    	SeqIO.write(record, corrected, 'fasta')

with open(specie2+".faa") as original, open(specie2+".faa2", 'w') as corrected:
    records = SeqIO.parse(specie2+".faa", 'fasta')
    for record in records:
    	record.description = record.id = record.id.split(':', 1)[0]
    	SeqIO.write(record, corrected, 'fasta')

with open(specie2+".fna") as original, open(specie2+".fna2", 'w') as corrected:
    records = SeqIO.parse(specie2+".fna", 'fasta')
    for record in records:
    	record.description = record.id = record.id.split(':', 1)[0]
    	SeqIO.write(record, corrected, 'fasta')

output_aa_sp1.close()
output_dna_sp1.close()
output_aa_sp2.close()
output_dna_sp2.close()

record_dict_sp1_aa = SeqIO.to_dict(SeqIO.parse(specie1+".faa2", "fasta"))
record_dict_sp2_aa = SeqIO.to_dict(SeqIO.parse(specie2+".faa2", "fasta"))
record_dict_sp1_dna = SeqIO.to_dict(SeqIO.parse(specie1+".fna2", "fasta"))
record_dict_sp2_dna = SeqIO.to_dict(SeqIO.parse(specie2+".fna2", "fasta"))

#Only keep sequences wich are present in the 2 species
with open(specie1+".faa", 'w') as corrected:
	for i in df2["Busco id"]:
		if i in record_dict_sp1_aa and record_dict_sp2_aa:
			SeqIO.write(record_dict_sp1_aa[i], corrected, 'fasta')

with open(specie1+".fna", 'w') as corrected:
	for i in df2["Busco id"]:
		if i in record_dict_sp1_aa:
			SeqIO.write(record_dict_sp1_dna[i], corrected, 'fasta')

with open(specie2+".faa", 'w') as corrected:
	for i in df1["Busco id"]:
		if i in record_dict_sp2_aa:
			SeqIO.write(record_dict_sp2_aa[i], corrected, 'fasta')

with open(specie2+".fna", 'w') as corrected:
	for i in df1["Busco id"]:
		if i in record_dict_sp2_dna:
			SeqIO.write(record_dict_sp2_dna[i], corrected, 'fasta')

os.remove(specie1+".faa2") 
os.remove(specie2+".faa2") 
os.remove(specie1+".fna2") 
os.remove(specie2+".fna2")

output_aa_sp1.close()
output_dna_sp1.close()
output_aa_sp2.close()
output_dna_sp2.close()


print("Number of complete Busco sequences in the specie",specie1,":",Number_complet_sp1)
print("Number of complete Busco sequences in the specie",specie2,":",Number_complet_sp2)
print("Number of complete Busco sequences present in both species: ",len(SeqIO.to_dict(SeqIO.parse(specie1+".faa", "fasta"))))
