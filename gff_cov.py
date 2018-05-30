import pandas as pd
import itertools
from itertools import *
from Bio import SeqIO
import numpy as np
import argparse
import argparse
import os, shutil
from shutil import copyfile
import fileinput
import glob

"""

Output files (11): If you only want to keep the filter p-value + cov, please uncomment the corresponding lines. 

-Augustus_clustering.tab : Outut table of Augustus predicted genes after passing through the clustering filter (max pident within each cluster)

Final output with candidates genes predicted by Augustus after passing through the filter p-value:
-candidates_genes_pvalue_sp1.tab (for specie 1) (optional)
-candidates_genes_pvalue_sp2.tab (for specie 2) (optional)
-candidates_genes_pvalue_s1_sp2.tab (for both species) (all informations)

Final output with candidates genes predicted by Augustus after passing through the filter p-value and cov.
-candidates_genes_pvalue_cov_0035.tab (optional)
-candidates_genes_pvalue_cov_0042.tab (optional)
-candidates_genes_pvalue_cov_s1_sp2.tab (for both species) (all informations)

Final fasta output with candidates genes predicted by Augustus after passing through the filter p-value (this file will be usefull for the blast against the nr db step).
-candidates_aa_pvalue_0035.fasta
-candidates_aa_pvalue_0042.fasta
-candidates_dna_pvalue_0035.fasta
-candidates_dna_pvalue_0042.fasta


"""

parser = argparse.ArgumentParser(description='Merge en keep clusters with the best pident')
parser.add_argument("-d1", "--dN_dS_Busco", help="introduce the dN_dS output of Busco genes")
parser.add_argument("-d2", "--dN_dS_Augustus", help="introduce the dN_dS output of Augustus genes")
parser.add_argument("-s1", "--specie1", help="introduce the name of the specie 1")
parser.add_argument("-s2", "--specie2", help="introduce the name of the specie 2")
parser.add_argument("-c1", "--cov_sp1", help="introduce the  Cov table of scaffold for the sp1")
parser.add_argument("-c2", "--cov_sp2", help="introduce the  Cov table of scaffold for the sp2")
parser.add_argument("-g1", "--gff_Augustus1", help="introduce the run augustus file of specie1 trained with sp1 (gff format)")
parser.add_argument("-g2", "--gff_Augustus2", help="introduce the run augustus file of specie1 trained with sp2 (gff format)")
parser.add_argument("-g3", "--gff_Augustus3", help="introduce the run augustus file of specie2 trained with sp2 (gff format)")
parser.add_argument("-g4", "--gff_Augustus4", help="introduce the run augustus file of specie2 trained with sp1 (gff format)")
parser.add_argument("-g5", "--gff_Busco1", help="introduce the run augustus of busco for specie1 (gff format)")
parser.add_argument("-g6", "--gff_Busco2", help="introduce the run augustus of busco for specie2 (gff format)")

args = parser.parse_args()

#Ex Usage: 
#python3 gff_cov.py -d1 dn_ds_Busco.out -d2 dn_ds_Augustus.out  -s1 0035 -s2 0042 -c1 cov_GC_0035.tab -c2 cov_GC_0042.tab -g1 run_augustus_0035_training_0035.out -g2 run_augustus_0035_training_0042.out -g3 run_augustus_0042_training_0042.out -g4 run_augustus_0042_training_0035.out -g5 gff_Busco_0035 -g6 gff_Busco_0042

# Variable that stores fasta sequences
dn_ds_Busco=args.dN_dS_Busco
dn_ds_Augustus=args.dN_dS_Augustus
specie1=args.specie1 #espèce 1
specie2=args.specie2 #espèce 2
gff_sp2_sp2=args.gff_Augustus3
gff_sp1_sp2=args.gff_Augustus2
gff_sp2_sp1=args.gff_Augustus4
gff_sp1_sp1=args.gff_Augustus1
cov_sp1=args.cov_sp1
cov_sp2=args.cov_sp2
gff_Busco_sp1=args.gff_Busco1
gff_Busco_sp2=args.gff_Busco2



"""
==============================================================================================================
=.     PREMIERE PARTIE
=
=. # Script qui prend en entrée la distribution des distances synonymes dN des gènes Busco (hypothèse nulle)
=. # Se base sur un seuil d'acceptation de transfert horizontal à 1% des dN des gènes Busco 
=. # Crée ensuite un fichier fasta acide aminé de ces séquences par espèce pour ensuite pouvoir les blaster 
=
==============================================================================================================
output_aa_file = open(specie1+'_aa.fasta','w')
output_aa_file2 = open(specie2+'_aa.fasta','w')

output_dna_file = open(specie1'+_dna.fasta','w')
output_dna_file2 = open(specie2+'_dna.fasta','w')

"""
#Count number of sequences present after clustering 
records = SeqIO.index("clusters1_aa.fasta", "fasta")
count_seq_after_cluster=len(records)


######################################
# I.     FILTERING  -- P-value.      #
######################################

#Open the Augustus predicted genes dataframe with distances informations 
dN_dS_not_filtred=pd.read_table(dn_ds_Augustus,header=0,sep="\t")

#Open the Busco dataframe with distances informations 
dS_busco=pd.read_table(dn_ds_Busco,header=0,sep=";")

dS_busco.columns = dS_busco.columns.str.replace('\s+', '_')  # in case there are multiple white spaces


dS_busco = dS_busco.drop(dS_busco[dS_busco.Mean_length < 750].index)

#p-value at 1%
p_value=dS_busco['dS'].quantile(q=0.01)

#Creat a subset containing candidates genes below p-value of 1%
dN_dS_not_filtred['dS'] < p_value
candidate_df=dN_dS_not_filtred[dN_dS_not_filtred['dS'] < p_value]


#candidate_df.to_csv("candidate_df",sep='\t')
#get the number of candidates genes under the p-value threshold
#print(candidate_df.shape)


#candidate_df=pd.read_table("candidate_df",header=0,sep="\t")

#Allow to transfomr the dN_dS file (one columns contains only sequences of sp1 and the other contains sequences of sp2), sorting GC and other informationq as well.
candidate_df.Length_seq_1, candidate_df.Length_seq_2 = np.where(candidate_df.seq1_id.str.contains('_'+specie1+'_'), candidate_df.Length_seq_1, candidate_df.Length_seq_2), np.where(candidate_df.seq1_id.str.contains('_'+specie2+'_'), candidate_df.Length_seq_1, candidate_df.Length_seq_2)
candidate_df.GC_content_seq1, candidate_df.GC_content_seq2 = np.where(candidate_df.seq1_id.str.contains('_'+specie1+'_'), candidate_df.GC_content_seq1, candidate_df.GC_content_seq2), np.where(candidate_df.seq1_id.str.contains('_'+specie2+'_'), candidate_df.GC_content_seq1, candidate_df.GC_content_seq2)
candidate_df.seq1_id, candidate_df.seq2_id = np.where(candidate_df.seq1_id.str.contains('_'+specie1+'_'), candidate_df.seq1_id, candidate_df.seq2_id), np.where(candidate_df.seq1_id.str.contains('_'+specie2+'_'), candidate_df.seq1_id, candidate_df.seq2_id)

#candidate_df.to_csv("dn_ds.out_test",sep='\t')

#Load the sequences comming from the cluster filtering and range them into ordered files per species

#candidate_df=pd.read_table("dn_ds.out_test",header=0,sep="\t")

#Creat a fasta file for each sequence species to perform a blast 
seq1_id=candidate_df["seq1_id"]
seq2_id=candidate_df["seq2_id"]

#Creat new fasta file with sequence passed through p_value filter:
output_aa_sp1 = open('candidates_aa_pvalue_'+specie1+'.fasta','w')
output_aa_sp2 = open('candidates_aa_pvalue_'+specie2+'.fasta','w')
output_dna_sp1 = open('candidates_dna_pvalue_'+specie1+'.fasta','w')
output_dna_sp2 = open('candidates_dna_pvalue_'+specie2+'.fasta','w')

#Open the fasta file with all sequences
record_dict_sp1_aa = SeqIO.to_dict(SeqIO.parse("clusters1_aa.fasta", "fasta"))
record_dict_sp2_aa = SeqIO.to_dict(SeqIO.parse("clusters2_aa.fasta", "fasta"))
record_dict_sp1_dna = SeqIO.to_dict(SeqIO.parse("clusters1_dna.fasta", "fasta"))
record_dict_sp2_dna = SeqIO.to_dict(SeqIO.parse("clusters2_dna.fasta", "fasta"))

#Amino Acide
for i in candidate_df["seq1_id"]:
    if i in record_dict_sp1_aa:
        SeqIO.write(record_dict_sp1_aa[i], output_aa_sp1, 'fasta')
    elif i in record_dict_sp2_aa:
        SeqIO.write(record_dict_sp2_aa[i], output_aa_sp1, 'fasta')

for i in candidate_df["seq2_id"]:
	if i in record_dict_sp1_aa:
		SeqIO.write(record_dict_sp1_aa[i], output_aa_sp2, 'fasta')
	elif i in record_dict_sp2_aa:
		SeqIO.write(record_dict_sp2_aa[i], output_aa_sp2, 'fasta')

#DNA
for i in candidate_df["seq1_id"]:
    if i in record_dict_sp1_dna:
        SeqIO.write(record_dict_sp1_dna[i], output_dna_sp1, 'fasta')
    elif i in record_dict_sp2_aa:
        SeqIO.write(record_dict_sp2_dna[i], output_dna_sp1, 'fasta')

for i in candidate_df["seq2_id"]:
	if i in record_dict_sp1_dna:
		SeqIO.write(record_dict_sp1_dna[i], output_dna_sp2, 'fasta')
	elif i in record_dict_sp2_aa:
		SeqIO.write(record_dict_sp2_dna[i], output_dna_sp2, 'fasta')

output_aa_sp1.close()
output_aa_sp2.close()
output_dna_sp1.close()
output_dna_sp2.close()

#Count how many sequences are filtred after p-value:
count_seq_after_p_value=0
for record in SeqIO.parse("candidates_aa_pvalue_"+specie2+".fasta", "fasta"):
    count_seq_after_p_value+=1

count_seq_after_p_value2=0
for record in SeqIO.parse("candidates_aa_pvalue_"+specie1+".fasta", "fasta"):
    count_seq_after_p_value2+=1


"""
===========================================================================
=
=.  DEUXIEME PARTIE
=
=.  Permet de récupérer le séquences après filtre de couverture et GC
=
=
=============================================================================
"""

####################################################################################################
#Getting augustus informations such as the scaffold number where are present the candidates genes  #
####################################################################################################

#sequences of specie sp2_sp2

liste=["scaf_name","source","feature","start","end","score","strand","frame","gene"]
gene_info=pd.read_csv(gff_sp2_sp2,comment='"',sep='\s+',header=None,names=liste)

scaf_info=pd.read_csv(cov_sp2,sep='\t')
scaf_info.scaf_name = scaf_info.scaf_name.str.replace(' ', '_')
ggf_sp2_sp2=pd.merge(gene_info, scaf_info, on='scaf_name')
ggf_sp2_sp2=ggf_sp2_sp2[ggf_sp2_sp2.feature == 'transcript']

#sequences of specie sp2_sp1
liste=["scaf_name","source","feature","start","end","score","strand","frame","gene"]
gene_info=pd.read_csv(gff_sp2_sp1,comment='"',sep='\s+',header=None,names=liste)

scaf_info=pd.read_csv(cov_sp2,sep='\t')
scaf_info.scaf_name = scaf_info.scaf_name.str.replace(' ', '_')
ggf_sp2_sp1=pd.merge(gene_info, scaf_info, on='scaf_name')
ggf_sp2_sp1=ggf_sp2_sp1[ggf_sp2_sp1.feature == 'transcript']

#sequences of specie sp1_sp1
liste=["scaf_name","source","feature","start","end","score","strand","frame","gene"]
gene_info=pd.read_csv(gff_sp1_sp1,comment='"',sep='\s+',header=None,names=liste)

scaf_info=pd.read_csv(cov_sp1,sep='\t')
scaf_info.scaf_name = scaf_info.scaf_name.str.replace(' ', '_')
ggf_sp1_sp1=pd.merge(gene_info, scaf_info, on='scaf_name')
ggf_sp1_sp1=ggf_sp1_sp1[ggf_sp1_sp1.feature == 'transcript']

#sequences of specie sp1_sp2
liste=["scaf_name","source","feature","start","end","score","strand","frame","gene"]
gene_info=pd.read_csv(gff_sp1_sp2,comment='"',sep='\s+',header=None,names=liste)

scaf_info=pd.read_csv(cov_sp1,sep='\t')
scaf_info.scaf_name = scaf_info.scaf_name.str.replace(' ', '_')
ggf_sp1_sp2=pd.merge(gene_info, scaf_info, on='scaf_name')
ggf_sp1_sp2=ggf_sp1_sp2[ggf_sp1_sp2.feature == 'transcript']

####################################################################################################
#Getting augustus informations such as the scaffold number where are present the candidates genes  #
####################################################################################################

#Busco sequences of specie sp1
liste=["scaf_name","source","feature","start","end","score","strand","frame","gene"]
gene_info_sp1=pd.read_csv(gff_Busco_sp1,comment='"',sep='\s+',header=None,names=liste)
scaf_info_sp1=pd.read_csv(cov_sp1,sep='\t')
scaf_info_sp1.scaf_name = scaf_info_sp1.scaf_name.str.replace(' ', '_')
gff_Busco_sp1=pd.merge(gene_info_sp1, scaf_info_sp1, on='scaf_name')
gff_Busco_sp1=gff_Busco_sp1[gff_Busco_sp1.feature == 'transcript']
gff_Busco_sp1= gff_Busco_sp1[["gene","scaf_name","start","end","cov_depth","GC"]]	
#gff_Busco_sp1.to_csv("Busco_"+specie1+"_cov_depth.txt",sep='\t')

#Busco sequences of specie sp2
liste=["scaf_name","source","feature","start","end","score","strand","frame","gene"]
gene_info_sp2=pd.read_csv(gff_Busco_sp2,comment='"',sep='\s+',header=None,names=liste)

scaf_info_sp2=pd.read_csv(cov_sp2,sep='\t')
scaf_info_sp2.scaf_name = scaf_info_sp2.scaf_name.str.replace(' ', '_')
gff_Busco_sp2=pd.merge(gene_info_sp2, scaf_info_sp2, on='scaf_name')
gff_Busco_sp2=gff_Busco_sp2[gff_Busco_sp2.feature == 'transcript']

gff_Busco_sp2=gff_Busco_sp2.drop(columns=["source","feature","score","strand","frame"])
gff_Busco_sp2= gff_Busco_sp2[["gene","scaf_name","start","end","cov_depth","GC"]]	
#gff_Busco_sp2.to_csv("Busco_"+specie2+"_cov_depth.txt",sep='\t')

#################################################################################################
#.  Getting the min and max distribution of busco "Gc and cov-depth" desired (5% quantile).     #
#################################################################################################

#If you want to add a GC filtering (but be carefull, we are looking for HGT, then the GC content of a recent transfered gene won't have the same GC as the host genome)
"""Busco_sp2_GC_max=gff_Busco_sp2['GC'].quantile(q=0.975)
#print("Busco_sp2_GC_max: ",Busco_sp2_GC_max)
Busco_sp2_GC_min=gff_Busco_sp2['GC'].quantile(q=0.025)
#print("Busco_sp2_GC_min: ",Busco_sp2_GC_min)
Busco_sp2_GC_mean=gff_Busco_sp2['GC'].mean()
#print("Busco_sp2_GC_mean: ",Busco_sp2_GC_mean)

Busco_sp1_GC_max=gff_Busco_sp1['GC'].quantile(q=0.975)
#print("Busco_sp1_GC_max: ",Busco_sp1_GC_max)
Busco_sp1_GC_min=gff_Busco_sp1['GC'].quantile(q=0.025)
#print("Busco_sp1_GC_min: ",Busco_sp1_GC_min)
Busco_sp1_GC_mean=gff_Busco_sp1['GC'].mean()
#print("Busco_sp1_GC_mean: ",Busco_sp1_GC_mean)"""


#Cov-depth filtering
Busco_sp2_cov_max=gff_Busco_sp2['cov_depth'].quantile(q=0.975)
#print("Busco_sp2_cov_max: ",Busco_sp2_cov_max)
Busco_sp2_cov_min=gff_Busco_sp2['cov_depth'].quantile(q=0.025)
#print("Busco_sp2_cov_min: ",Busco_sp2_cov_min)
Busco_sp2_cov_mean=gff_Busco_sp2['cov_depth'].mean()
#print("Busco_sp2_cov_mean: ",Busco_sp2_cov_mean)

Busco_sp1_cov_max=gff_Busco_sp1['cov_depth'].quantile(q=0.975)
#print("Busco_sp1_cov_max: ",Busco_sp1_cov_max)
Busco_sp1_cov_min=gff_Busco_sp1['cov_depth'].quantile(q=0.025)
#print("Busco_sp1_cov_min: ",Busco_sp1_cov_min)
Busco_sp1_cov_mean=gff_Busco_sp1['cov_depth'].mean()
#print("Busco_sp1_cov_mean: ",Busco_sp1_cov_mean)

#The following part is necessary only if you have renamed you augusuts predicted genes into their Busco names
"""
#Correct busco name to augustus name in candidates transfered genes
sp1_aa="candidates_aa_sp1.fasta"
sp1_aa_corrected="candidates_aa_sp1_corrected.fasta"
sp2_aa="candidates_aa_sp2.fasta"
sp2_aa_corrected="candidates_aa_sp2_corrected.fasta"



with open(sp2_aa) as original, open(sp2_aa_corrected, 'w') as corrected:
    records = SeqIO.parse(sp2_aa, 'fasta')
    for record in records:
        if record.id == "EOG090X0IUZ_sp1_sp1_1":
            print(record.id)
            record.id = "g7044.t1_sp1_sp1"
            record.description = "g7044.t1_sp1_sp1"

        SeqIO.write(record, corrected, 'fasta')
  

with open(sp1_aa) as original, open(sp1_aa_corrected, 'w') as corrected:
    records = SeqIO.parse(sp1_aa, 'fasta')
    for record in records:
        if record.id == "EOG090X0IUZ_sp2_sp2_1":
            print(record.id)
            record.id = "g5713.t1_sp2_sp2"
            record.description = "g5713.t1_sp2_sp2"

        if record.id == "EOG090X01YQ_sp2_sp2_2":
            print(record.id)
            record.id = "g13545.t1_sp2_sp2"
            record.description = "g13545.t1_sp2_sp2" # <- Add this line

        SeqIO.write(record, corrected, 'fasta')


"""
#Next part allow to creat 2 dataframes with informations such : [gene,scaf_name,start,end,cov_depth,GC], one df for each specie c"andidates_genes_sp1" and "candidates_genes_sp2".
# Then it also concatenate those two dataframe into one "candidates_genes"

df2_sp1 = pd.DataFrame(columns=("scaf_name","source","feature","start","end","score","strand","frame","gene"))
df2_sp2 = pd.DataFrame(columns=("scaf_name","source","feature","start","end","score","strand","frame","gene"))

#####adding to dataframe informations about cov percentage

for record in SeqIO.parse("candidates_aa_pvalue_"+specie1+".fasta", 'fasta'): #pars the dataframe
	gene_name=str(record.id.split('_', maxsplit=1)[0])
	if(record.id[record.id.index("_"):]) == "_"+specie1+"_"+specie1: #if the number_number_ = _sp1_sp1, then,
		df2_sp1=df2_sp1.append(ggf_sp1_sp1[ggf_sp1_sp1['gene']== gene_name])
		df2_sp1.loc[df2_sp1['gene'] == gene_name, 'gene'] += "_"+specie1+"_"+specie1

	if(record.id[record.id.index("_"):]) == "_"+specie1+"_"+specie2: #if the number_number_ = _sp1_sp1, then,
		df2_sp1=df2_sp1.append(ggf_sp1_sp2[ggf_sp1_sp2['gene']== gene_name])
		df2_sp1.loc[df2_sp1['gene'] == gene_name, 'gene'] += "_"+specie1+"_"+specie2

	if(record.id[record.id.index("_"):]) == "_"+specie2+"_"+specie1: #if the number_number_ = _sp1_sp1, then,
		df2_sp2=df2_sp2.append(ggf_sp2_sp1[ggf_sp2_sp1['gene']== gene_name])
		df2_sp2.loc[df2_sp2['gene'] == gene_name, 'gene'] += "_"+specie2+"_"+specie1

	if(record.id[record.id.index("_"):]) == "_"+specie2+"_"+specie2: #if the number_number_ = _sp1_sp1, then,
		df2_sp2=df2_sp2.append(ggf_sp2_sp2[ggf_sp2_sp2['gene']== gene_name])
		df2_sp2.loc[df2_sp2['gene'] == gene_name, 'gene'] += "_"+specie2+"_"+specie2


#####adding to dataframe informations about cov percentage

for record in SeqIO.parse("candidates_aa_pvalue_"+specie2+".fasta", 'fasta'): #pars the dataframe
	gene_name=str(record.id.split('_', maxsplit=1)[0])
	if(record.id[record.id.index("_"):]) == "_"+specie1+"_"+specie1: #if the number_number_ = _sp1_sp1, then,
		df2_sp1=df2_sp1.append(ggf_sp1_sp1[ggf_sp1_sp1['gene']== gene_name])
		df2_sp1.loc[df2_sp1['gene'] == gene_name, 'gene'] += "_"+specie1+"_"+specie1

	if(record.id[record.id.index("_"):]) == "_"+specie1+"_"+specie2: #if the number_number_ = _sp1_sp1, then,
		df2_sp1=df2_sp1.append(ggf_sp1_sp2[ggf_sp1_sp2['gene']== gene_name])
		df2_sp1.loc[df2_sp1['gene'] == gene_name, 'gene'] += "_"+specie1+"_"+specie2

	if(record.id[record.id.index("_"):]) == "_"+specie2+"_"+specie1: #if the number_number_ = _sp1_sp1, then,
		df2_sp2=df2_sp2.append(ggf_sp2_sp1[ggf_sp2_sp1['gene']== gene_name])
		df2_sp2.loc[df2_sp2['gene'] == gene_name, 'gene'] += "_"+specie2+"_"+specie1

	if(record.id[record.id.index("_"):]) == "_"+specie2+"_"+specie2: #if the number_number_ = _sp1_sp1, then,
		df2_sp2=df2_sp2.append(ggf_sp2_sp2[ggf_sp2_sp2['gene']== gene_name])
		df2_sp2.loc[df2_sp2['gene'] == gene_name, 'gene'] += "_"+specie2+"_"+specie2


df2_sp2=df2_sp2.drop(columns=["source","feature","score","strand","frame"])
df2_sp2= df2_sp2[["gene","scaf_name","start","end","cov_depth","GC"]]	
df2_sp1=df2_sp1.drop(columns=["source","feature","score","strand","frame"])
df2_sp1= df2_sp1[["gene","scaf_name","start","end","cov_depth","GC"]]	
df2_sp2.to_csv("candidates_genes_pvalue_"+specie2+".tab",sep='\t')
df2_sp1.to_csv("candidates_genes_pvalue_"+specie1+".tab",sep='\t')

#Merge the two dataframe
df2=df2_sp1.append(df2_sp2)



#Creat a new dataframe with all informations of these paired sequences afetr a p-value filter:
augustus_cluster=pd.read_table("Augustus_clustering.tab",header=0,sep='\t')
record_dict_sp1 = SeqIO.to_dict(SeqIO.parse('candidates_aa_pvalue_'+specie1+'.fasta', "fasta"))

dataframe = augustus_cluster[augustus_cluster['seq1_id'].isin(record_dict_sp1) | augustus_cluster['seq2_id'].isin(record_dict_sp1)]

cov_tab_sp1=pd.read_table("candidates_genes_pvalue_"+specie1+".tab",header=0,sep='\t')

cov_tab_sp2=pd.read_table("candidates_genes_pvalue_"+specie2+".tab",header=0,sep='\t')

data = dataframe.merge(cov_tab_sp1, left_on='seq1_id', right_on='gene')
data = data.merge(cov_tab_sp2, left_on='seq2_id', right_on='gene')
data=data.drop(columns=["Unnamed: 0_x", "Unnamed: 0.1", "Unnamed: 0_y", "gene_x", 'Unnamed: 0','gene_x','start_y','end_x','start_x','end_y','gene_y','GC_x'])

data.columns = ['seq1_id','seq2_id','dN','dS','Dist_third_pos','Dist_brute',
 'Length_seq_1','Length_seq_2','GC_content_seq1','GC_content_seq2',
 'Mean_length','scaf_name_seq1','cov_depth_seq1','GC_scaff_seq1','scaf_name_seq2','cov_depth_seq2','GC_scaff_seq2']

data= data[['seq1_id','seq2_id','dN','dS','Dist_third_pos','Dist_brute','Mean_length',
 'Length_seq_1','Length_seq_2','GC_content_seq1','GC_content_seq2','scaf_name_seq1','GC_scaff_seq1','cov_depth_seq1','scaf_name_seq2','GC_scaff_seq2','cov_depth_seq2']]
data.to_csv("candidates_genes_pvalue_"+specie1+"_"+specie2+".tab",sep='\t')



######################################
# II     FILTERING  -- COV -- GC.    #
######################################

#Filtering on specie 2
candidates_sp2=pd.read_csv("candidates_genes_pvalue_"+specie2+".tab",sep='\t')
candidates_sp2 = candidates_sp2[candidates_sp2.cov_depth < Busco_sp2_cov_max]
candidates_sp2 = candidates_sp2[candidates_sp2.cov_depth > Busco_sp2_cov_min]
#candidates_sp2 = candidates_sp2[candidates_sp2.GC < Busco_sp2_GC_max]
#candidates_sp2 = candidates_sp2[candidates_sp2.GC > Busco_sp2_GC_min]
candidates_sp2.to_csv("candidates_genes_pvalue_cov_"+specie2+".tab",sep='\t')

#Fitering on specie 1
candidates_sp1=pd.read_csv("candidates_genes_pvalue_"+specie1+".tab",sep='\t')
candidates_sp1 = candidates_sp1[candidates_sp1.cov_depth < Busco_sp1_cov_max]
candidates_sp1 = candidates_sp1[candidates_sp1.cov_depth > Busco_sp1_cov_min]
#candidates_sp1 = candidates_sp1[candidates_sp1.GC < Busco_sp1_GC_max]
#candidates_sp1 = candidates_sp1[candidates_sp1.GC > Busco_sp1_GC_min]
candidates_sp1.to_csv("candidates_genes_pvalue_cov_"+specie1+".tab",sep='\t')



#################################
#. Sorte the dN_dS file         #
#################################
#Open the dataframe with distances informations 
dN_dS=pd.read_table(dn_ds_Augustus,header=0,sep="\t")

#Allow to transfomr the dN_dS file (one columns contains only sequences of sp1 and the other contains sequences of sp2), sorting GC and other informationq as well.
dN_dS.Length_seq_1, dN_dS.Length_seq_2 = np.where(dN_dS.seq1_id.str.contains('_'+specie1+'_'), dN_dS.Length_seq_1, dN_dS.Length_seq_2), np.where(dN_dS.seq1_id.str.contains('_'+specie2+'_'), dN_dS.Length_seq_1, dN_dS.Length_seq_2)
dN_dS.GC_content_seq1, dN_dS.GC_content_seq2 = np.where(dN_dS.seq1_id.str.contains('_'+specie1+'_'), dN_dS.GC_content_seq1, dN_dS.GC_content_seq2), np.where(dN_dS.seq1_id.str.contains('_'+specie2+'_'), dN_dS.GC_content_seq1, dN_dS.GC_content_seq2)
dN_dS.seq1_id, dN_dS.seq2_id = np.where(dN_dS.seq1_id.str.contains('_'+specie1+'_'), dN_dS.seq1_id, dN_dS.seq2_id), np.where(dN_dS.seq1_id.str.contains('_'+specie2+'_'), dN_dS.seq1_id, dN_dS.seq2_id)

dN_dS.to_csv("Augustus_clustering.tab",sep='\t')

#_________________________

#Because in a pairs, one seq can passe the filter while the other does not, it is important to keep only paired genes wich are keeping on the two previous filtering:

gene_name_sp1=[]
for i in candidates_sp1["gene"]:
	gene_name_sp1.append(i)

gene_name_sp2=[]
for i in candidates_sp2["gene"]:
	gene_name_sp2.append(i)

df4 = pd.DataFrame(columns=dN_dS.columns)
for index, row in dN_dS.iterrows():

    if row['seq1_id'] in gene_name_sp1 and row['seq2_id'] in gene_name_sp2:
        df4 = df4.append(row, ignore_index=True)

candidates_sp1=[]
for i in df4["seq1_id"]:
    candidates_sp1.append(i)

candidates_sp1 = pd.DataFrame({'seq1_id':candidates_sp1})
#candidates_sp1.to_csv("conserved_candidates_"+specie1,sep='\t')

candidates_sp2=[]
for i in df4["seq2_id"]:
    candidates_sp2.append(i)

candidates_sp2 = pd.DataFrame({'seq2_id':candidates_sp2})
#candidates_sp2.to_csv("conserved_candidates_"+specie2,sep='\t')

#Get output after filtering

Dataframe_HGT_candidates_filtered=pd.merge(dN_dS, df4, how='inner', on=['seq1_id', 'seq2_id'])

Dataframe_HGT_candidates_filtered.to_csv("candidates_genes_pvalue_cov_"+specie1+"_"+specie2+".tab",sep='\t')

#If you want to keep these files, uncomment this code
os.remove("candidates_genes_pvalue_"+specie2+".tab") 
os.remove("candidates_genes_pvalue_"+specie1+".tab") 
os.remove("candidates_genes_pvalue_cov_"+specie2+".tab") 
os.remove("candidates_genes_pvalue_cov_"+specie1+".tab")



print(" =====================================")
print(" =. Criteria retained : ","\n","=====================================","\n",
	"Max coverage",specie2,": ",round(Busco_sp2_cov_max,3),"\n","Min coverage",specie2,": ",round(Busco_sp2_cov_min,3),"\n","Max coverage",specie1,": "
	,round(Busco_sp1_cov_max,3),"\n","Min coverage",specie1,": ",round(Busco_sp1_cov_min,3),"\n")

print("Number of sequence pairs after clustering silix: ", count_seq_after_cluster)
print("Number of sequence pairs after p-value threshold of", p_value," :", count_seq_after_p_value ,"(Sequences that will be analyzed against the nr database)")
print("Number of sequences after purification 'cover':","\n", "-sp_",specie2,":", candidates_sp2.shape[0],"\n","-sp_",specie1,":",candidates_sp1.shape[0],"\n","-total:", candidates_sp2.shape[0]+candidates_sp1.shape[0])
print("Number of pairs of sequences after purification 'cover':", Dataframe_HGT_candidates_filtered.shape[0])






