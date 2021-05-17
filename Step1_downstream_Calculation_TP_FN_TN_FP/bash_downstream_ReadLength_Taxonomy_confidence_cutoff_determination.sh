#!/bin/bash

WorkDir=/mnt/home1/qinglong/Test_SimulatedAmplicon_NCBI16SRefSeq_ReadLength_taxonomy/ReadLength_taxa_accuracy/Step1_downstream_Calculation_TP_FN_TN_FP

cd $WorkDir

for file in $WorkDir/*_*.txt  #V1V3_forward.txt prepared from Windows
do
	FileName=`basename $file .txt`
	#extract every 4 columns
	cut -d$'\t' -f 1-4 $file > ''$FileName'-100nt.txt'
	cut -d$'\t' -f 5-8 $file > ''$FileName'-150nt.txt'
	cut -d$'\t' -f 9-12 $file > ''$FileName'-170nt.txt'
	cut -d$'\t' -f 13-16 $file > ''$FileName'-200nt.txt'
	cut -d$'\t' -f 17-20 $file > ''$FileName'-250nt.txt'
	cut -d$'\t' -f 21-24 $file > ''$FileName'-300nt.txt'
	cut -d$'\t' -f 25-28 $file > ''$FileName'-350nt.txt'
	cut -d$'\t' -f 29-32 $file > ''$FileName'-400nt.txt'
	cut -d$'\t' -f 33-36 $file > ''$FileName'-450nt.txt'
	cut -d$'\t' -f 37-40 $file > ''$FileName'-Original.txt'  #Tab inconsistent between Linux and Windows, so add anything at 41th column will sove the proble for calculation at species level 

	#calculate accuracy at genus level
	for SubFile in $WorkDir/*_*-*.txt
	do
		SubFileName=`basename $SubFile .txt`
		for threshold in 70 80 85 90 95 100
		do
			#note the default OFS of awk is space-delimited
			TP=`cat $SubFile | cut -d$'\t' -f 1,3 | awk -F '\t' '$2=="Match" {print $1,$2}' | awk '{if ($1 >= '$threshold') print $1,$2}' | wc -l`
			FN=`cat $SubFile | cut -d$'\t' -f 1,3 | awk -F '\t' '$2=="Match" {print $1,$2}' | awk '{if ($1 < '$threshold') print $1,$2}' | wc -l`
			TN=`cat $SubFile | cut -d$'\t' -f 1,3 | awk -F '\t' '$2=="Not-Match" {print $1,$2}' | awk '{if ($1 < '$threshold') print $1,$2}' | wc -l`
			FP=`cat $SubFile | cut -d$'\t' -f 1,3 | awk -F '\t' '$2=="Not-Match" {print $1,$2}' | awk '{if ($1 >= '$threshold') print $1,$2}' | wc -l`
			echo -e "$SubFileName\t$threshold\t$TP\t$FN\t$TN\t$FP"
		done
	done > temp1
	echo $'Region-Orientation-Length\tThreshold\tTP\tFN\tTN\tFP' | cat - temp1 > ''$FileName'.taxonomic.accuracy.genus.txt'
	rm temp1

	#calculate accuracy at species level
        for SubFile in $WorkDir/*_*-*.txt
        do
                SubFileName=`basename $SubFile .txt`
                for threshold in 30 40 50 60 70 80
                do
                        #note the default OFS of awk is space-delimited
                        TP=`cat $SubFile | cut -d$'\t' -f 2,4 | awk -F '\t' '$2=="Match" {print $1,$2}' | awk '{if ($1 >= '$threshold') print $1,$2}' | wc -l`
                        FN=`cat $SubFile | cut -d$'\t' -f 2,4 | awk -F '\t' '$2=="Match" {print $1,$2}' | awk '{if ($1 < '$threshold') print $1,$2}' | wc -l`
                        TN=`cat $SubFile | cut -d$'\t' -f 2,4 | awk -F '\t' '$2=="Not-Match" {print $1,$2}' | awk '{if ($1 < '$threshold') print $1,$2}' | wc -l`
                        FP=`cat $SubFile | cut -d$'\t' -f 2,4 | awk -F '\t' '$2=="Not-Match" {print $1,$2}' | awk '{if ($1 >= '$threshold') print $1,$2}' | wc -l`
                        echo -e "$SubFileName\t$threshold\t$TP\t$FN\t$TN\t$FP"
                done
        done > temp2
        echo $'Region-Orientation-Length\tThreshold\tTP\tFN\tTN\tFP' | cat - temp2 > ''$FileName'.taxonomic.accuracy.species.txt'
        rm temp2

	rm *_*-*.txt
done
