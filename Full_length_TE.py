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
from TEinsertion_bamConverter import bamConverter
import numpy as np
import warnings
warnings.filterwarnings('ignore')


pd.set_option("display.max_column",40)

parser=argparse.ArgumentParser()
parser.add_argument("-f1","--file1")
parser.add_argument("-f2","--file2")
parser.add_argument("-pre","--Prefix")
args=parser.parse_args()


TEfile=args.file1
Gefile=args.file2
preFix=args.Prefix

def getMultiFragment(f):
	f=f.sort_values(["Readname","ReadStart_chr","ReadEnd_chr"])
	f2_1=f.drop_duplicates(["Readname","TE","ReadStart_TE","ReadStart_TE"],keep="first")
	f2_2=f.drop_duplicates(["Readname","TE","ReadStart_TE","ReadStart_TE"],keep="last")
	new_f2=f2_1.merge(f2_2,on=["Readname","TE","ReadStart_TE","ReadEnd_TE"],how="inner")
	new_f2=new_f2.drop(["TE_Start_y","TE_End_y","ReadLen_y","Strand_TE_y","ReadLen_y","flanking_x","flanking_y"],axis=1)
	new_f2.columns=["Readname","ReadLen","TE","TE_Start","TE_End","ReadStart_TE","ReadEnd_TE","Strand_TE","chr_x","chr_Start_x","chr_End_x","ReadStart_chr_x","ReadEnd_chr_x","Strand_chr_x","dis1_x","dis2_x","chr_y","chr_Start_y", "chr_End_y","ReadStart_chr_y","ReadEnd_chr_y","Strand_chr_y","dis1_y","dis2_y"]
	new_f2["dis1"]=abs(new_f2["ReadEnd_chr_x"]-new_f2["ReadStart_TE"])
	new_f2["dis2"]=abs(new_f2["ReadStart_chr_y"]-new_f2["ReadEnd_TE"])
	new_f2=new_f2.loc[(new_f2["dis1"]<=1000) & (new_f2["dis2"]<=1000)]

	return new_f2
	

def fullLength(TEfile,Gefile):
	TE=pd.read_table(TEfile)
	Gefile=pd.read_table(Gefile)
	full=TE.loc[(TE["RefStart"]==1) & (TE["RefEnd"]==7469)]
	G=Gefile.loc[Gefile["Readname"].isin(list(full["Readname"]))]
	allAlig=full.merge(G,on="Readname",how="inner")
	allAlig.columns=["TE","TE_Start","TE_End","Readname","ReadLen","ReadStart_TE","ReadEnd_TE","Strand_TE","chr","chr_Start","chr_End","ReadLen_y","ReadStart_chr","ReadEnd_chr","Strand_chr"]
	allAlig=allAlig[["Readname","ReadLen","TE","TE_Start","TE_End","ReadStart_TE","ReadEnd_TE","Strand_TE","chr","chr_Start","chr_End","ReadStart_chr","ReadEnd_chr","Strand_chr"]]
	return allAlig

def NonSplitCopies(allAlignment):
	allAlignment["flanking"]=False
	mask1=(allAlignment["ReadStart_chr"]<=100) & (allAlignment["ReadEnd_chr"]>=allAlignment["ReadLen"]-100)
	mask2=(allAlignment["ReadStart_chr"]-allAlignment["ReadStart_TE"]<=-100) & (allAlignment["ReadEnd_chr"]-allAlignment["ReadEnd_TE"]>=100)	
	allAlignment.loc[mask1,"flanking"]=True
	allAlignment.loc[mask2,"flanking"]=True
	RefIns=allAlignment.loc[allAlignment["flanking"]==True]
	RefIns["t1"]=allAlignment["ReadStart_chr"]-allAlignment["ReadStart_TE"]
	RefIns["t2"]=allAlignment["ReadEnd_chr"]-allAlignment["ReadEnd_TE"]

	RefIns["J1"]=0
	RefIns.loc[RefIns["Strand_chr"]=="+","J1"]=RefIns["chr_Start"]+(RefIns["ReadStart_TE"]-RefIns["ReadStart_chr"])
	RefIns.loc[RefIns["Strand_chr"]=="-","J1"]=RefIns["chr_Start"]+(RefIns["ReadEnd_chr"]-RefIns["ReadEnd_TE"])
	RefIns["J2"]=0
	RefIns.loc[RefIns["Strand_chr"]=="+","J2"]=RefIns["chr_End"]-(RefIns["ReadEnd_chr"]-RefIns["ReadEnd_TE"])
	RefIns.loc[RefIns["Strand_chr"]=="-","J2"]=RefIns["chr_End"]-(RefIns["ReadStart_TE"]-RefIns["ReadStart_chr"])
	RefIns["J"]=RefIns["J2"]-RefIns["J1"]
	RefIns=RefIns[["Readname","ReadLen","TE","TE_Start","TE_End","ReadStart_TE","ReadEnd_TE","Strand_TE","chr","chr_Start","chr_End","ReadStart_chr","ReadEnd_chr","Strand_chr","J1","J2"]]
	RefIns["type"]="nonSplit"
	return RefIns

def SplitCopies(allAlign,NonSplit):
	SplCop=allAlign.loc[~allAlign["Readname"].isin(list(NonSplit["Readname"]))].sort_values(["Readname"])
	SplCop["flanking"]=False
	mask1=(SplCop["ReadStart_chr"]<=100) | (SplCop["ReadEnd_chr"]>=SplCop["ReadLen"]-100)
	mask2=(SplCop["ReadStart_chr"]-SplCop["ReadStart_TE"]<=-100) | (SplCop["ReadEnd_chr"]-SplCop["ReadEnd_TE"]>=100)
	SplCop.loc[mask1,"flanking"]=True	
	SplCop.loc[mask2,"flanking"]=True
	SplCop=SplCop.loc[SplCop["flanking"]==True]
	SplCop["dis1"]=abs(SplCop["ReadEnd_chr"]-SplCop["ReadStart_TE"])
	SplCop["dis2"]=abs(SplCop["ReadEnd_TE"]-SplCop["ReadStart_chr"])
	SplCop=SplCop.loc[(SplCop["dis1"]<=1000) | (SplCop["dis2"]<=1000)]
	s=SplCop.groupby(["Readname"]).filter(lambda x: len(x)==1)
	s=s[["Readname","ReadLen","TE","TE_Start","TE_End","ReadStart_TE","ReadEnd_TE","Strand_TE","chr","chr_Start","chr_End","ReadStart_chr","ReadEnd_chr","Strand_chr"]]
	s["J1"]=0
	s.loc[(s["ReadStart_TE"]<s["ReadStart_chr"])& (s["Strand_chr"]=="+"),"J1"]=s["chr_Start"]
	s.loc[(s["ReadStart_TE"]<s["ReadStart_chr"])& (s["Strand_chr"]=="-"),"J1"]=s["chr_End"]
	s.loc[(s["ReadStart_TE"]>s["ReadStart_chr"])& (s["Strand_chr"]=="+"),"J1"]=s["chr_End"]
	s.loc[(s["ReadStart_TE"]>s["ReadStart_chr"])& (s["Strand_chr"]=="-"),"J1"]=s["chr_Start"]
	s["J2"]=s["J1"]
	s["type"]="OneSide"
	
	single=SplCop.loc[~SplCop["Readname"].isin(list(s["Readname"]))]
	new_f2=pd.DataFrame(columns=["Readname","ReadLen","TE","TE_Start","TE_End","ReadStart_TE","ReadEnd_TE","Strand_TE","chr_x","chr_Start_x","chr_End_x","ReadStart_chr_x","ReadEnd_chr_x","Strand_chr_x","dis1_x","dis2_x","chr_y","chr_Start_y", "chr_End_y","ReadStart_chr_y","ReadEnd_chr_y","Strand_chr_y","dis1_y","dis2_y","dis1","dis2"])
	for i in list(set(single["Readname"])):
		sub_f=single.loc[single["Readname"]==i]
		I=sub_f.index.values.tolist()
		i_1=0
		while i_1<len(I)-1:
			i_2=i_1+1
			while i_2 <len(I):
				sub_sub_f=sub_f.loc[[I[i_1],I[i_2]]]
				sub_sub_f_res=getMultiFragment(sub_sub_f)
				if sub_sub_f_res.shape[0]>0:
					new_f2=new_f2.append(sub_sub_f_res)
					i_2=i_2+1
				i_2=i_2+1
			i_1=i_1+1
	new_f2=new_f2.drop(["dis1_x","dis2_x","dis1_y","dis2_y","dis1","dis2"],axis=1)
	return s,new_f2

def GetJunction(NonSplit,s,Split):
	df1=NonSplit[["Readname","chr","J1","J2","type"]]

	c2=Split.loc[Split["chr_x"]==Split["chr_y"]].sort_values(["Readname"])
	#c2=c2.groupby(["Readname","ReadStart_TE","ReadEnd_TE"]).filter(lambda x: len(x)==1)
	c2=c2.reset_index()
	c2["J1"]=0
	c2.loc[c2["Strand_chr_x"]=="+","J1"]=c2["chr_End_x"]
	c2.loc[c2["Strand_chr_x"]=="-","J1"]=c2["chr_Start_x"]
	c2["J2"]=0
	c2.loc[c2["Strand_chr_y"]=="+","J2"]=c2["chr_Start_y"]
	c2.loc[c2["Strand_chr_y"]=="-","J2"]=c2["chr_End_y"]
	c2["type"]="split"
	c2["chr"]=c2["chr_x"]
	c2=c2.sort_values(["Readname","ReadStart_chr_x","ReadEnd_chr_x","ReadStart_chr_y","ReadEnd_chr_y"])
	c2["JunDis"]=abs(c2["J1"]-c2["J2"])
	c2=c2.sort_values(["Readname","JunDis"])
	c2=c2.drop_duplicates(["Readname","ReadStart_TE","ReadEnd_TE"],keep="first")
	
	df2=c2[["Readname","chr","J1","J2","type"]]
	
	df3=s[["Readname","chr","J1","J2","type"]]

	df=df1.append(df2)
	df=df.append(df3)
	df=df.reset_index()
	df["S"]=0
	df["E"]=0
	df.loc[df["J1"]<=df["J2"],"S"]=df["J1"]
	df.loc[df["J1"]<=df["J2"],"E"]=df["J2"]
	df.loc[df["J1"]>df["J2"],"S"]=df["J2"]
	df.loc[df["J1"]>df["J2"],"E"]=df["J1"]
	
	df=df[["chr","S","E","J1","J2","Readname","type"]]
	df=df.sort_values(["chr","S","E"])
	df.to_csv(preFix+"_mapped_copies.tsv",header=None,index=None,sep="\t")
	bedtools="bedtools cluster -d 1000 -i %s > %s"%(preFix+"_mapped_copies.tsv",preFix+"_mapped_copies.tsv.cluster.bed")
	os.system(bedtools)
	os.system("rm %s"%(preFix+"_mapped_copies.tsv"))
	bedout=pd.read_table(preFix+"_mapped_copies.tsv.cluster.bed",header=None)	
	most_coor={}
	for clu in set(bedout[7]):
		sub_clu=bedout.loc[bedout[7]==clu]
		coor=sub_clu[1].value_counts().idxmax()
		most_coor[clu]=coor
	bedout[8]=bedout[7].apply(lambda x: most_coor[x])
	
	bedout.columns=["chr","S","E","J1","J2","Readname","type","clusterID","cluster_coor"]
	bedout["cluster_coor2"]=bedout["cluster_coor"]+1
	bedout=bedout[["chr","cluster_coor","cluster_coor2","J1","J2","Readname","type","clusterID"]]
	bedout=bedout.loc[bedout["Readname"].isnull()==False]
	bedout.to_csv(preFix+"_mapped_copies.tsv.cluster.bed",index=None,sep="\t")
	return bedout

def MultiFrag(Split,bedout):
	c3=Split.loc[Split["chr_x"]!=Split["chr_y"]].sort_values(["Readname"])
	multi=c3.loc[~c3["Readname"].isin(list(bedout["Readname"]))]
	multi["dis1"]=abs(multi["ReadEnd_chr_x"]-multi["ReadStart_TE"])
	multi["dis2"]=abs(multi["ReadStart_chr_y"]-multi["ReadEnd_TE"])
	multi=multi.sort_values(["Readname","dis1","dis2"])
	multi=multi.drop_duplicates(["Readname","ReadStart_TE","ReadEnd_TE"],keep="first")
	multi=multi.reset_index()
	chrList=["chr2R","chr2L","chr3R","chr3L","chr4","chrX","chrY"]
	multi=multi.loc[(multi["chr_x"].isin(chrList)) & (multi["chr_y"].isin(chrList))]
	multi["J1"]=0
	multi.loc[multi["Strand_chr_x"]=="+","J1"]=multi["chr_End_x"]
	multi.loc[multi["Strand_chr_x"]=="-","J1"]=multi["chr_Start_x"]
	multi["J2"]=0
	multi.loc[multi["Strand_chr_y"]=="+","J2"]=multi["chr_Start_y"]
	multi.loc[multi["Strand_chr_y"]=="-","J2"]=multi["chr_End_y"]
	multi["J1.1"]=multi["J1"]+1
	#print(multi.shape)
	#print(multi.drop_duplicates(["Readname"],keep="first").shape)
	#print(multi)

	links=multi[["chr_x","J1","chr_y","J2","Readname"]]	
	links.to_csv(preFix+"_links.tsv",header=None,index=None,sep="\t")
	df1=multi[["chr_x","J1","J1.1","Readname"]]
	df1=df1.sort_values(["chr_x","J1"])
	df1.to_csv(preFix+"_Unmapped_copies.tsv",header=None,index=None,sep="\t")
	bedtools="bedtools cluster -d 1000 -i %s > %s"%(preFix+"_Unmapped_copies.tsv",preFix+"_Unmapped_copies.tsv.cluster.bed")
	os.system(bedtools)
	os.system("rm %s"%(preFix+"_Unmapped_copies.tsv"))
	f=pd.read_table(preFix+"_Unmapped_copies.tsv.cluster.bed",header=None)
	#print(preFix,multi.drop_duplicates(["Readname"],keep="first").shape[0],max(f[4]))
	l=list(multi["Readname"])+list(bedout["Readname"])
	others=c3.loc[~c3["Readname"].isin(l)]
	#print(others)
	#print(others.drop_duplicates(["Readname"],keep="first").shape)
	return multi

allAlign=fullLength(TEfile,Gefile)
#print(allAlign.drop_duplicates(["Readname"],keep="first").shape)
NonSplit=NonSplitCopies(allAlign)
#print(NonSplit.drop_duplicates(["Readname"],keep="first").shape)
#print(NonSplit[0:10])
s,SplitCop=SplitCopies(allAlign,NonSplit)

bedout=GetJunction(NonSplit,s,SplitCop)
c=bedout.drop_duplicates(["Readname"],keep="first").groupby(["clusterID"],as_index=False).count()[["clusterID","Readname"]].sort_values(["Readname"],ascending=[False])
#print(c.shape,c["Readname"].sum())
c_m=c.loc[c["Readname"]>=2]
c_u=c.loc[c["Readname"]==1]
#print(c_m.shape)
#print(c_u.shape)

c=bedout.drop_duplicates(["Readname"],keep="first")
c_m_file=c.loc[c["clusterID"].isin(c_m["clusterID"])].drop_duplicates(["clusterID"],keep="first")[["chr","cluster_coor","clusterID"]]
c_m_file.to_csv(preFix+"_mapped_multi_loc.tsv",header=None,sep="\t",index=None)

c_u_file=c.loc[c["clusterID"].isin(c_u["clusterID"])].drop_duplicates(["clusterID"],keep="first")[["chr","cluster_coor","clusterID"]]
c_u_file.to_csv(preFix+"_mapped_single_loc.tsv",header=None,sep="\t",index=None)


#GetMul(SplitCop)
#print(allAlign)
multi=MultiFrag(SplitCop,bedout)


#print(multi)
#print(multi.shape)

total=allAlign.drop_duplicates(["Readname"],keep="first").shape[0]
assign=bedout.drop_duplicates(["Readname"],keep="first").shape[0]
location=bedout.drop_duplicates(["clusterID"],keep="first").shape[0]
print(preFix,str(total),str(assign),str(location))


def seperate(Pre):
	f=pd.read_table(Pre+"_Unmapped_copies.tsv.cluster.bed",header=None)
	c=f.drop_duplicates([3],keep="first").groupby([4],as_index=False).count()[[4,3]]
	c_m=c.loc[c[3]>=2]
	c_u=c.loc[c[3]==1]
#	print(c_m.shape)
#	print(c_u.shape)
	c_m_file=f.loc[f[4].isin(c_m[4])]
	c_u_file=f.loc[f[4].isin(c_u[4])]
	l=pd.read_table(Pre+"_links.tsv",header=None,sep="\t")
	l_m=l.loc[l[4].isin(list(c_m_file[3]))].sort_values([0,1,2,3])
	l_u=l.loc[l[4].isin(list(c_u_file[3]))].sort_values([0,1,2,3])
#	print(l_m.shape)
#	print(l_u.shape)
	l_m.to_csv(Pre+"_links.multi.tsv",index=None,header=None,sep="\t")
	l_u.to_csv(Pre+"_links.single.tsv",index=None,header=None,sep="\t")

#seperate("barcode21")
#seperate("barcode22")

