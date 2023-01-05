import pandas as pd
import os
import os.path
import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
from cigar import Cigar
from collections import Counter
import pysam
import warnings
warnings.filterwarnings('ignore')


pd.set_option("display.max_column",40)

parser=argparse.ArgumentParser()
parser.add_argument("-Ta","--TE_bam")
parser.add_argument("-Ga","--Genome_bam")
parser.add_argument("-pName","--Prefix")
parser.add_argument("-flex","--flexibility",default=100)

args=parser.parse_args()

Ta=args.TE_bam
Ga=args.Genome_bam
pName=args.Prefix
fl=args.flexibility


def combine(TEfile,Genomefile):
	TE=pd.read_table(TEfile,header=None,sep=" ")
	Genome=pd.read_table(Genomefile,header=None,sep=" ")
	TE=TE[range(0,9)]
	Genome=Genome[range(0,9)]
	TE_columns=["Readname","ReadLen","ReadStart_TE","ReadEnd_TE","Strand_TE","TE_Name","TElen","TEstart","TEend"]
	Genome_columns=["Readname","ReadLen","ReadStart_REF","ReadEnd_REF","Strand_REF","REF_Name","REFlen","REFstart","REFend"]
	TE.columns=TE_columns
	Genome.columns=Genome_columns
	print(TE.shape)
	print(TE[0:10])
	print(Genome.shape)
	print(Genome[0:10])
	combined=pd.merge(TE,Genome,how="inner",on=["Readname","ReadLen"])

	combined["overlap"]=False
	# remove full comverage
	full_cover=((combined["ReadStart_REF"]<combined["ReadStart_TE"]) & (combined["ReadEnd_REF"]>combined["ReadStart_TE"]))
	combined.loc[full_cover,"overlap"]=True
	remove=combined.loc[combined["overlap"]==True]
	combined=combined.loc[~(combined["Readname"].isin(list(remove["Readname"])))]

	print(combined.shape)
	print(combined[0:10])
	print(combined.drop_duplicates("Readname",keep="first").shape)
	
	combined["overlap"]=True
	conf=((combined["ReadStart_REF"]>=combined["ReadEnd_TE"]-fl) | (combined["ReadEnd_REF"]<=combined["ReadStart_TE"]+fl))
	combined.loc[conf,"overlap"]=False
	
	f=combined.loc[combined["overlap"]==False]
	f["dis1"]=abs(f["ReadEnd_REF"]-f["ReadStart_TE"])
	f["dis2"]=abs(f["ReadStart_REF"]-f["ReadEnd_TE"])
	f=f.loc[(f["dis1"]<=100) | (f["dis2"]<=100)]
	print(f.shape)
	print(f[0:10])
	print(f.drop_duplicates("Readname",keep="first").shape)
	f=f.drop(["overlap"],axis=1)
	f.to_csv(pName+"_filtered.tsv",index=None,sep="\t")
#	
#combine(Ta,Ga)

def single(filteredFile):
	f=pd.read_table(filteredFile,header=0,sep="\t")
	f=f.drop_duplicates(["Readname","ReadStart_REF","ReadEnd_REF"])
	single=f.groupby(["Readname","TE_Name"],as_index=False).filter(lambda x: len(x)==1)
	single["Junc_1"]=0
	single.loc[(single["Strand_REF"]=="-")&(single["dis1"]<=100),"Junc_1"]=single["REFstart"].apply(int)
	single.loc[(single["Strand_REF"]=="-")&(single["dis2"]<=100),"Junc_1"]=single["REFend"].apply(int)
	single.loc[(single["Strand_REF"]=="+")&(single["dis1"]<=100),"Junc_1"]=single["REFend"].apply(int)
	single.loc[(single["Strand_REF"]=="+")&(single["dis2"]<=100),"Junc_1"]=single["REFstart"].apply(int)
	
	print(single[0:10])
	print(single.shape)
	print(single.drop_duplicates(["Readname"],keep="first").shape)

#single(pName+"_filtered.tsv")

def double(filteredFile):
	f=pd.read_table(filteredFile,header=0,sep="\t")
	f=f.drop_duplicates(["Readname","ReadStart_REF","ReadEnd_REF"])
	double=f.groupby(["Readname","TE_Name"],as_index=False).filter(lambda x: len(x)==2)

	double["Junc_1"]=0
	double["Junc_2"]=0
	double.loc[(double["Strand_REF"]=="-")&(double["dis1"]<=100),"Junc_1"]=double["REFstart"].apply(int)
	double.loc[(double["Strand_REF"]=="+")&(double["dis1"]<=100),"Junc_1"]=double["REFend"].apply(int)
	double.loc[(double["Strand_REF"]=="+")&(double["dis2"]<=100),"Junc_2"]=double["REFstart"].apply(int)
	double.loc[(double["Strand_REF"]=="-")&(double["dis2"]<=100),"Junc_2"]=double["REFend"].apply(int)
	print(double[0:10])
	print(double.shape)
	print(double.drop_duplicates(["Readname"],keep="first").shape)
	

double(pName+"_filtered.tsv")


#def getMultiFragment(f):
#  f2_1=f.drop_duplicates(["Readname","TEname","ReadStart_TE","ReadStart_TE"],keep="first")
#  f2_2=f.drop_duplicates(["Readname","TEname","ReadStart_TE","ReadStart_TE"],keep="last")
#  new_f2=f2_1.merge(f2_2,on=["Readname","TEname","ReadStart_TE","ReadEnd_TE"],how="inner")
#  new_f2=new_f2.drop(["ReadLen_y","Strand_TE_y","TEstart_y","TEend_y","TElength_y"],axis=1)
#  new_f2=new_f2[["Readname","ReadLen_x","ReadStart_TE","ReadEnd_TE","Strand_TE_x","TEname","TElength_x","TEstart_x","TEend_x","RefName_x","RefStart_x","RefEnd_x","ReadStart_ref_x","ReadEnd_ref_x","Strand_ref_x","RefName_y","RefStart_y","RefEnd_y","ReadStart_ref_y","ReadEnd_ref_y","Strand_ref_y"]]
#
#  new_f2.columns=["Readname","ReadLen","ReadStart_TE","ReadEnd_TE","Strand_TE","TEname","TElength_x","TEstart_x","TEend_x","RefName_x","RefStart_x","RefEnd_x","ReadStart_ref_x","ReadEnd_ref_x","Strand_ref_x","RefName_y","RefStart_y","RefEnd_y","ReadStart_ref_y","ReadEnd_ref_y","Strand_ref_y"]
#  
#  new_f2=new_f2.loc[new_f2["RefName_x"]==new_f2["RefName_y"]]
#  new_f2["dis1"]=abs(new_f2["ReadEnd_ref_x"]-new_f2["ReadStart_TE"])
#  new_f2["dis2"]=abs(new_f2["ReadStart_ref_y"]-new_f2["ReadEnd_TE"])
#  new_f2=new_f2.loc[(new_f2["dis1"]<=100) & (new_f2["dis2"]<=100)]
#  new_f2=new_f2.loc[(new_f2["ReadStart_ref_x"]<=fl) & (new_f2["ReadEnd_ref_y"]>=new_f2["ReadLen"]-fl)]
#  return new_f2
#
def get_insertions(filteredFile):
	f=pd.read_table(filteredFile,header=0,sep="\t")
	f=f.sort_values(["Readname","TEname","ReadStart_TE","ReadEnd_TE","RefName","ReadStart_ref","ReadEnd_ref"])
	f_d=f.drop_duplicates(["Readname","TEname","ReadStart_TE","ReadEnd_TE","ReadStart_ref","ReadEnd_ref"],keep="first")
	f1=f_d.groupby(["Readname","TEname","ReadStart_TE","ReadEnd_TE"]).filter(lambda x: len(x)==1)
	f2=f_d.groupby(["Readname","TEname","ReadStart_TE","ReadEnd_TE"]).filter(lambda x: len(x)>=2)
	new_f2=pd.DataFrame(columns=["Readname","ReadLen","ReadStart_TE","ReadEnd_TE","Strand_TE","TEname","TElength_x","TEstart_x","TEend_x","RefName_x","RefStart_x","RefEnd_x","ReadStart_ref_x","ReadEnd_ref_x","Strand_ref_x","RefName_y","RefStart_y","RefEnd_y","ReadStart_ref_y","ReadEnd_ref_y","Strand_ref_y"])
	for i in list(set(f2["Readname"])):
		sub_f=f.loc[f["Readname"]==i]
		I=sub_f.index.values.tolist()
		i_1=0
		while i_1<len(I)-1:
			i_2=i_1+1
			while i_2 <len(I):
				sub_sub_f=sub_f.loc[[I[i_1],I[i_2]]]
				sub_sub_f_res=getMultiFragment(sub_sub_f)
				if sub_sub_f_res.shape[0]>0:
					new_f2.append(sub_sub_f_res)
					i_2=i_2+1
				i_2=i_2+1
			i_1=i_1+1
	f1.to_csv(pName+"_org_single.tsv",index=None,sep="\t")
	new_f2.to_csv(pName+"_org_double.tsv",index=None,sep="\t")

	
def re_countDouble(double_df):
	f=pd.read_table(double_df,header=0,sep="\t")
	f["type"]="double"
	if f.shape[0]==0:
		return f
	else:
		f["Junction_1"]=0
		f["Junction_2"]=0
		f.loc[f["Strand_ref_x"]=="+","Junction_1"]=f["RefEnd_x"]
		f.loc[f["Strand_ref_x"]=="-","Junction_1"]=f["RefStart_x"]
		f.loc[f["Strand_ref_y"]=="+","Junction_2"]=f["RefStart_y"]
		f.loc[f["Strand_ref_y"]=="-","Junction_2"]=f["RefEnd_y"]
		f=f.loc[(f["Junction_2"]-f["Junction_1"]<=100) & (f["Junction_2"]-f["Junction_1"]>=-100)]
		return f

def re_countSingle(orgFile):
	f=pd.read_table(orgFile,header=0,sep="\t")
	f["sum"]=f["ReadEnd_ref"]-f["ReadStart_ref"]+1+f["ReadEnd_TE"]-f["ReadStart_TE"]+1
	f=f.loc[f["sum"]/f["ReadLen"]>=0.9]
	f["dis1"]=abs(f["ReadEnd_ref"]-f["ReadStart_TE"])
	f["dis2"]=abs(f["ReadStart_ref"]-f["ReadEnd_TE"])
	f["junction_1"]=0
	f["junction_2"]=0
	f=f.loc[(f["dis1"]<=100) | (f["dis2"]<=100)]
	f.loc[(f["Strand_ref"]=="-")&(f["dis1"]<=100),"junction_1"]=f["RefStart"]
	f.loc[(f["Strand_ref"]=="-")&(f["dis2"]<=100),"junction_1"]=f["RefEnd"]
	f.loc[(f["Strand_ref"]=="+")&(f["dis1"]<=100),"junction_1"]=f["RefEnd"]
	f.loc[(f["Strand_ref"]=="+")&(f["dis2"]<=100),"junction_1"]=f["RefStart"]
#	
	f["min"]=0
	f["max"]=0
	f.loc[f["dis1"]<=100,"min"]=f["ReadStart_ref"]
	f.loc[f["dis1"]<=100,"max"]=f["ReadEnd_TE"]
	f.loc[f["dis2"]<100,"min"]=f["ReadStart_TE"]
	f.loc[f["dis2"]<100,"max"]=f["ReadEnd_ref"]
	f["TEdir"]=False
	f.loc[(f["min"]==f["ReadStart_TE"]) & (f["Strand_TE"]=="+")& (f["TEend"]>=f["TElength"]-100),"TEdir"]=True
	f.loc[(f["min"]==f["ReadStart_TE"]) & (f["Strand_TE"]=="-")& (f["TEstart"]<=100),"TEdir"]=True
	f.loc[(f["min"]==f["ReadStart_ref"])& (f["Strand_TE"]=="+")&(f["TEstart"]<=100),"TEdir"]=True
	f.loc[(f["min"]==f["ReadStart_ref"])& (f["Strand_TE"]=="-")&(f["TEend"]>=f["TElength"]-100),"TEdir"]=True
	f=f.loc[f["TEdir"]==True]
	f=f.loc[(f["min"]<=100) & (f["max"]>=f["ReadLen"]-100)]
	f=f.drop_duplicates(["Readname","ReadStart_TE","ReadEnd_TE"],keep="first")
	for i in ["RefName_y","RefStart_y","RefEnd_y","RefLen_y","ReadStart_ref_y","ReadEnd_ref_y","Strand_ref_y"]:
		f[i]="NA"
	
	f=f[["Readname","ReadLen","ReadStart_TE","ReadEnd_TE","Strand_TE","TEname","TElength","TEstart","TEend","RefName","RefStart","RefEnd","ReadStart_ref","ReadEnd_ref","Strand_ref","RefName_y","RefStart_y","RefEnd_y","ReadStart_ref_y","ReadEnd_ref_y","Strand_ref_y","dis1","dis2","junction_1","junction_2"]]
	f.columns=["Readname","ReadLen","ReadStart_TE","ReadEnd_TE","Strand_TE","TEname","TElength_x","TEstart_x","TEend_x","RefName_x","RefStart_x","RefEnd_x","ReadStart_ref_x","ReadEnd_ref_x","Strand_ref_x","RefName_y","RefStart_y","RefEnd_y","ReadStart_ref_y","ReadEnd_ref_y","Strand_ref_y","dis1","dis2","Junction_1","Junction_2"]
	f["type"]="single"
	return f
#


#AddTElength(Ta+"_AligTable.tsv")
#combine(Ta+"_AligTable.tsv_TE.tsv",Ga+"_AligTable.tsv")
#get_insertions(pName+"_filtered.tsv")
#
#f2=re_countDouble(pName+"_org_double.tsv")
#f1=re_countSingle(pName+"_org_single.tsv")
##
#if f2.shape[0]>0:
#	final=f1.append(f2)
#else:
#	final=f1
#
##
#final["cluster"]=final["Junction_1"].apply(lambda x: int(round(x,-2)))
#final_insertion=final.groupby(["RefName_x","cluster"],as_index=False).count()
#final_insertion=final.groupby(["RefName_x","Junction_1"],as_index=False).count()
#final_insertion["i"]=final_insertion["RefName_x"]+"_"+final_insertion["Junction_1"].apply(str)
#
#d=dict(zip(final_insertion["i"],final_insertion["Readname"]))
##
#final=final.sort_values(["RefName_x","Junction_1"])
#final["i"]=final["RefName_x"]+"_"+final["Junction_1"].apply(str)
#final["count"]=final["i"].apply(lambda x: d[x])
#final=final.drop(["i"],axis=1)
#final.to_csv(pName+"_Final_insertion.tsv",index=None,sep="\t")
#
#print(final.shape)
#print(final_insertion.shape)
