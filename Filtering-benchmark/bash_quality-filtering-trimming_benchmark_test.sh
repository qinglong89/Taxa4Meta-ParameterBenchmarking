#!/bin/bash

WorkDir=/mnt/home1/qinglong/Test_SimulatedAmplicon_NCBI16SRefSeq_ReadLength_taxonomy/Filtering-benchmark
THREADS=5

#MinLen=300  #for 454 reads
MinLen=200  #for Illumina reads

cd $WorkDir

#for InputSeq in $WorkDir/PRJNA48333_V1V3_all.fastq $WorkDir/PRJNA48333_V3V5_all.fastq
#for InputSeq in $WorkDir/UK_Illumina_V4_forward_all.fastq $WorkDir/UK_Illumina_V1V3_forward_all.fastq
for InputSeq in $WorkDir/UK_Illumina_V4_reverse_all.fastq $WorkDir/UK_Illumina_V1V3_reverse_all.fastq
do
	cd $WorkDir

	SeqName=`basename $InputSeq .fastq`

	mkdir $SeqName
	cd $SeqName

	#Average Q filtering
	cat $InputSeq | NanoFilt --quality 20 --length $MinLen > ''$SeqName'_AverageQ-20-'$MinLen'nt.fastq' & #keep in background

	#Minimum Q filtering
	/home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --fastq_filter $InputSeq \
								  --fastqout ''$SeqName'_MinimumQ-20-'$MinLen'nt.fastq' \
								  --fastq_minlen $MinLen \
								  --fastq_truncqual 20 \
								  --threads $THREADS & #keep in background

	#Maxmimum EE filtering
	/home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --fastq_filter $InputSeq \
        	                                                  --fastqout ''$SeqName'_MaxEE-1-'$MinLen'nt.fastq' \
                	                                          --fastq_minlen $MinLen \
                        	                                  --fastq_truncee 1 \
                                	                          --threads $THREADS & #keep in background

	#Sliding window filtering
	java -jar /home/qinglong/softwares/trimmomatic-binary-0.36/trimmomatic-0.36.jar SE -phred33 -threads $THREADS \
        	                                                                        $InputSeq \
                	                                                                ''$SeqName'_Slide-4-15-'$MinLen'nt.fastq' \
	                       	                                                        SLIDINGWINDOW:4:15 \
                                	                                                MINLEN:$MinLen & #keep in background


	java -jar /home/qinglong/softwares/trimmomatic-binary-0.36/trimmomatic-0.36.jar SE -phred33 -threads $THREADS \
											$InputSeq \
											''$SeqName'_Slide-20-25-'$MinLen'nt.fastq' \
											SLIDINGWINDOW:20:25 \
											MINLEN:$MinLen & #keep in background

	wait

	#Sliding window filtering + Maxmimum EE filtering
	/home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --fastq_filter ''$SeqName'_Slide-4-15-'$MinLen'nt.fastq' \
        	                                                  --fastqout ''$SeqName'_Slide-4-15-'$MinLen'nt_MaxEE-1.fastq' \
                	                                          --fastq_maxee 1 \
                        	                                  --threads $THREADS & #keep in background

	/home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --fastq_filter ''$SeqName'_Slide-20-25-'$MinLen'nt.fastq' \
        	                                                  --fastqout ''$SeqName'_Slide-20-25-'$MinLen'nt_MaxEE-1.fastq' \
                	                                          --fastq_maxee 1 \
                        	                                  --threads $THREADS & #keep in background
	wait


	#calculate the mean and stdev of the read length of the quality-passed reads
	for TrimmedFile in $PWD/*.fastq
	do
		TrimmedFileName=`basename $TrimmedFile .fastq`

		cat $TrimmedFile | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > ''$TrimmedFileName'_read_length_summary.txt'
		cat $TrimmedFile | awk '{if(NR%4==2) print length($1)}' | sort -n > ''$TrimmedFileName'_read_length.txt'

		AVERAGE=`awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' $TrimmedFile`
		STDEV=`awk '{sum+=$1; sumsq+=$1*$1}END{print sqrt(sumsq/NR - (sum/NR)**2)}' ''$TrimmedFileName'_read_length.txt'`
		TotalReads=`echo $(cat $TrimmedFile | wc -l)/4|bc`

		echo -e "$TrimmedFileName\t$AVERAGE\t$STDEV\t$TotalReads"
	done > ../''$SeqName'_stats_each_filtering.txt'

	#generate states for the sequences
	for TrimmedFile in $PWD/*.fastq
        do
		TrimmedFileName=`basename $TrimmedFile .fastq`

                /home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --fastq_eestats $TrimmedFile \
                                                                          --output ''$TrimmedFileName'_eestats.txt' \
                                                                          --threads $THREADS  & #keep in background
	done

	wait #until stats generation finished
done
