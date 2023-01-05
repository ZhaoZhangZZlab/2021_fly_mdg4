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
	f=f.sort_values(["Readname","TE_Name","ReadStart_REF","ReadEnd_REF"])
	f=f.drop_duplicates(["Readname","ReadStart_REF"],keep="last")
	f=f.drop_duplicates(["Readname","ReadEnd_REF"],keep="first")
	single=f.groupby(["Readname","TE_Name"],as_index=False).filter(lambda x: len(x)==1)
	single["Junc_1"]=-1
	single["Junc_2"]=-1
	single.loc[(single["Strand_REF"]=="-")&(single["dis1"]<=100),"Junc_1"]=single["REFstart"].apply(int)
	single.loc[(single["Strand_REF"]=="+")&(single["dis1"]<=100),"Junc_1"]=single["REFend"].apply(int)
	single.loc[(single["Strand_REF"]=="-")&(single["dis2"]<=100),"Junc_2"]=single["REFend"].apply(int)
	single.loc[(single["Strand_REF"]=="+")&(single["dis2"]<=100),"Junc_2"]=single["REFstart"].apply(int)
	
	return single

single_df=single(pName+"_filtered.tsv")

def double(filteredFile):
	f=pd.read_table(filteredFile,header=0,sep="\t")
	f=f.sort_values(["Readname","TE_Name","ReadStart_REF","ReadEnd_REF"])
	f=f.drop_duplicates(["Readname","ReadStart_REF"],keep="last")
	f=f.drop_duplicates(["Readname","ReadEnd_REF"],keep="first")
	double=f.groupby(["Readname","TE_Name"],as_index=False).filter(lambda x: len(x)==2)

	double["Junc_1"]=-1
	double["Junc_2"]=-1
	double.loc[(double["Strand_REF"]=="-")&(double["dis1"]<=100),"Junc_1"]=double["REFstart"].apply(int)
	double.loc[(double["Strand_REF"]=="+")&(double["dis1"]<=100),"Junc_1"]=double["REFend"].apply(int)
	double.loc[(double["Strand_REF"]=="+")&(double["dis2"]<=100),"Junc_2"]=double["REFstart"].apply(int)
	double.loc[(double["Strand_REF"]=="-")&(double["dis2"]<=100),"Junc_2"]=double["REFend"].apply(int)
	return double	

double_df=double(pName+"_filtered.tsv")


def single_reads(single_df):
	single_df.loc[single_df["dis1"]<=100,"left_refName"]=single_df["REF_Name"]
	single_df.loc[single_df["dis1"]<=100,"left_refStart"]=single_df["REFstart"].apply(int)
	single_df.loc[single_df["dis1"]<=100,"left_refEnd"]=single_df["REFend"].apply(int)
	single_df.loc[single_df["dis1"]<=100,"right_refName"]=-1
	single_df.loc[single_df["dis1"]<=100,"right_refStart"]=-1
	single_df.loc[single_df["dis1"]<=100,"right_refEnd"]=-1
	single_df.loc[single_df["dis1"]<=100,"Read_leftStart"]=single_df["ReadStart_REF"].apply(int)
	single_df.loc[single_df["dis1"]<=100,"Read_leftEnd"]=single_df["ReadEnd_REF"].apply(int)
	single_df.loc[single_df["dis1"]<=100,"Read_leftStrand"]=single_df["Strand_REF"]
	single_df.loc[single_df["dis1"]<=100,"Read_rightStart"]=-1
	single_df.loc[single_df["dis1"]<=100,"Read_rightEnd"]=-1
	single_df.loc[single_df["dis1"]<=100,"Read_rightStrand"]=-1

	single_df.loc[single_df["dis2"]<=100,"left_refName"]=-1
	single_df.loc[single_df["dis2"]<=100,"left_refStart"]=-1
	single_df.loc[single_df["dis2"]<=100,"left_refEnd"]=-1
	single_df.loc[single_df["dis2"]<=100,"right_refName"]=single_df["REF_Name"]
	single_df.loc[single_df["dis2"]<=100,"right_refStart"]=single_df["REFstart"].apply(int)
	single_df.loc[single_df["dis2"]<=100,"right_refEnd"]=single_df["REFend"].apply(int)
	single_df.loc[single_df["dis2"]<=100,"Read_leftStart"]=-1
	single_df.loc[single_df["dis2"]<=100,"Read_leftEnd"]=-1
	single_df.loc[single_df["dis2"]<=100,"Read_leftStrand"]=-1
	single_df.loc[single_df["dis2"]<=100,"Read_rightStart"]=single_df["ReadStart_REF"].apply(int)
	single_df.loc[single_df["dis2"]<=100,"Read_rightEnd"]=single_df["ReadEnd_REF"].apply(int)
	single_df.loc[single_df["dis2"]<=100,"Read_rightStrand"]=single_df["Strand_REF"]


	single_out=single_df[["Readname","ReadLen","TE_Name","TElen","left_refName","left_refStart","left_refEnd","TEstart","TEend","right_refName","right_refStart","right_refEnd","Read_leftStart","Read_leftEnd","Read_leftStrand","ReadStart_TE","ReadEnd_TE","Strand_TE","Read_rightStart","Read_rightEnd","Read_rightStrand","Junc_1","Junc_2"]]
	
	return single_out
single_out=single_reads(single_df)


def double_reads(double_df):
	double_df=double_df.sort_values(["Readname","TE_Name","ReadStart_REF","ReadEnd_REF"])
	double_df_1=double_df.drop_duplicates(["Readname","TE_Name","ReadStart_TE","ReadEnd_TE"],keep="first")
	double_df_2=double_df.drop_duplicates(["Readname","TE_Name","ReadStart_TE","ReadEnd_TE"],keep="last")

	double_df_1_out=single_reads(double_df_1)
	double_df_2_out=single_reads(double_df_2)
	double_out=pd.merge(double_df_1,double_df_2,on=["Readname","ReadLen","TE_Name","TElen","ReadStart_TE","ReadEnd_TE","Strand_TE","TEstart","TEend"],how="inner")
	remove1=double_out["Strand_REF_x"]!=double_out["Strand_REF_y"]
	remove2=double_out["REF_Name_x"]!=double_out["REF_Name_y"]
	remove=double_out.loc[remove1 | remove2]
	double_out=double_out.loc[~double_out["Readname"].isin(list(remove["Readname"]))]
	
	double_out=double_out[["Readname","ReadLen","TE_Name","TElen","REF_Name_x","REFstart_x","REFend_x","TEstart","TEend","REF_Name_y","REFstart_y","REFend_y","ReadStart_REF_x","ReadStart_REF_x","Strand_REF_x","ReadStart_TE","ReadEnd_TE","Strand_TE","ReadStart_REF_y","ReadEnd_REF_y","Strand_REF_y","Junc_1_x","Junc_2_y"]]
	
	columns=["Readname","ReadLen","TE_Name","TElen","left_refName","left_refStart","left_refEnd","TEstart","TEend","right_refName","right_refStart","right_refEnd","Read_leftStart","Read_leftEnd","Read_leftStrand","ReadStart_TE","ReadEnd_TE","Strand_TE","Read_rightStart","Read_rightEnd","Read_rightStrand","Junc_1","Junc_2"]
	double_out.columns=columns
	return double_out

double_out=double_reads(double_df)

def appendResult(single_out,double_out):
	single_out["conf"]="single"
	double_out["conf"]="double"
	
	f=single_out.append(double_out)
	full_TE=(f["TEstart"]<=fl) & (f["TEend"]>=f["TElen"]-fl)
	truncated=((f["TEstart"]<=fl) & (f["TEend"]<=f["TElen"]-fl)) | ((f["TEstart"]>=fl) & (f["TEend"]>=f["TElen"]-fl))
	frag=(f["TEstart"]>fl) & (f["TEend"]<f["TElen"]-fl)
	f.loc[full_TE,"TEconf"]="FL_TE"
	f.loc[truncated,"TEconf"]="trunc_TE"
	f.loc[frag,"TEconf"]="frag_TE"
	f.to_csv(pName+"_InsReads.tsv",index=None,sep="\t")
	return f
outf=appendResult(single_out,double_out)

def genomeLocation(outfile):
	f=pd.read_table(outfile)
	f["coor1"]=-1
	f["coor2"]=-1
	f.loc[f["Junc_2"]==-1,"coor1"]=f["Junc_1"].apply(lambda x: round(x,-2))
	f.loc[f["Junc_1"]==-1,"coor1"]=f["Junc_2"].apply(lambda x: round(x,-2))
	f.loc[(f["Junc_1"]!=-1)&(f["Junc_2"]!=-1) & (f["Junc_1"]<=f["Junc_2"]), "coor1"]=f["Junc_1"].apply(lambda x: round(x,-2))
	f.loc[(f["Junc_1"]!=-1)&(f["Junc_2"]!=-1) & (f["Junc_1"]>f["Junc_2"]), "coor1"]=f["Junc_2"].apply(lambda x: round(x,-2))
	f.loc[(f["Junc_1"]!=-1)&(f["Junc_2"]!=-1) & (f["Junc_1"]<=f["Junc_2"]), "coor2"]=f["Junc_2"].apply(lambda x: round(x,-2))
	f.loc[(f["Junc_1"]!=-1)& (f["Junc_2"]!=-1)& (f["Junc_1"]>=f["Junc_2"]), "coor2"]=f["Junc_1"].apply(lambda x: round(x,-2))
	f.loc[f["left_refName"]!="-1","refName"]=f["left_refName"]
	f.loc[f["left_refName"]=="-1","refName"]=f["right_refName"]
	f=f.loc[f["TEconf"]=="FL_TE"]
	f=f.groupby(["refName","coor1"],as_index=False).count()
	print(f.shape)
	print(f[0:10])


genomeLocation(pName+"_InsReads.tsv")
