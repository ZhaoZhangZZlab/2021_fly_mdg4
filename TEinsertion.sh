#!/bin/bash

read=$1
readName=$(basename $read);
NanoStat --fastq $read -n $readName".stats"
porechop -i $read --extra_end_trim 0 --discard_middle -o $readName".pre.fastq"
ref=$2
TE=$3
refName=$(basename $ref);
TEname=$(basename $TE)
minimap2 -x map-ont $ref $reads -Y -t 16 | awk '{for(i=1;i<=9;i++) printf $i" "; print ""}' FS='\t' > $readName"_"$refName".paf";
minimap2 -x map-ont $TE $reads -Y -t 16 | awk '{for(i=1;i<=9;i++) printf $i" "; print ""}' FS='\t' > $readName"_"$TEName".paf";
python3 /data/zhanglab/Weijia_Su/Git/TEinsertion/TEinsertion_20230102.py -Ga $readName"_"$refName".paf" -Ta $readName"_"$TEName".paf" -pName $readName 

