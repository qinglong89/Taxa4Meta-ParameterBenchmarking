#!/bin/bash

WorkDir=/mnt/home1/qinglong/Test_SimulatedAmplicon_NCBI16SRefSeq_ReadLength_taxonomy/Clustering_Denoising_accuracy
THREADS=30

for AmpliconDir in $WorkDir/NCBI_16S-rRNA-RefSeq_*_AbundanceSimulated_*_amplicon
do
	AmpliconDirName=`basename $AmpliconDir`
	cd $AmpliconDir

#	mkdir DADA2_v1.8

	#generate random quality score for each nucleotide of each sequence from a range (here, from 30 to 42)
	#DADA2 version 1.8 used ShortRead Package internal for parsing and processing fastq file
	#My test: for faking the quality score of simulated data, DADA2 v1.8 only accept quality scores of ASCII_BASE=64 style
	#Just need to change the fastq format in "fasta2fastq-better.py" which adopted Bio.Seq.IO for parsing fasta and fastq
#	for fastafile in $PWD/NCBI_16S-rRNA-RefSeq_*_AbundanceSimulated_*_*.fasta
#	do
#		fastafilename=`basename $fastafile .fasta`
#		python /mnt/home1/qinglong/Test_SimulatedAmplicon_NCBI16SRefSeq_ReadLength_taxonomy/fasta2fastq-better.py \
#			$fastafile ./DADA2_v1.8/''$fastafilename'.fastq' &	#keep in background
#	done
#	wait

	#above step generate fastq files to be used by DADA2
	#run the DADA2 commands for each directory, make sure you have R libraries installed as indicated in the R script
	cd DADA2_v1.8
#	Rscript --vanilla /mnt/home1/qinglong/Test_SimulatedAmplicon_NCBI16SRefSeq_ReadLength_taxonomy/bash_DADA2_script.R $AmpliconDir/DADA2_v1.8

        #parse the dada2 feature tables
#        cat DADA2_ASV-feature-table.txt | datamash transpose --no-strict > DADA2_ASV-feature-table_transposed.txt

	#parse the dada2 ASVs sequences: remove first line, convert tab file to fasta file, remove " character
#	source /home/DataAnalysis/miniconda2/bin/activate base  #activate conda base environment
#	cat DADA2_ASV-seq.txt | sed 1d | seqkit tab2fx | seqkit seq -w 0 | sed 's/"//g' > DADA2_ASV-seq_fixed.fasta
#	source /home/DataAnalysis/miniconda2/bin/deactivate base  #deactivate conda base environment

        #assign taxonomy by BLCA with NCBI 16S Microbial DB (downloaded on October 6, 2018)
	mkdir BLCA_taxonomy

        #BLCA current version just use single thread for the script
        #split input into several new files of relatively even size; then run the script on each individual file,
        #your will need to have several bash script for each subseq and later concatenate outputs together
        source /home/DataAnalysis/miniconda2/bin/activate base  #activate conda base environment
        seqkit split --by-part $THREADS DADA2_ASV-seq_fixed.fasta --out-dir Splits   #split sequences into N parts for downstream parallel processing
        source /home/DataAnalysis/miniconda2/bin/deactivate base

        for SplitSeq in $PWD/Splits/*.fasta
        do
                SplitSeqName=`basename $SplitSeq .fasta`
                python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i $SplitSeq \
									    --cvrset 0.99 --iset 99 --proc 1 \
                                                                            -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
                                                                            -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
                                                                            -o $PWD/BLCA_taxonomy/''$SplitSeqName'_BLCA_out.txt' &  #keep in background
        done

        wait  #wait the completion of all the background runs

        cat ./BLCA_taxonomy/*_BLCA_out.txt > ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_temp.txt
        rm ./BLCA_taxonomy/*_BLCA_out.txt  #remove the BLCA results for each sequence split

        rm -rf Splits  #remove the sequence splits

        ###unfornatunately, BLCA sometimes fails to annotate some sequences from the input during large run
        cut -d$'\t' -f 1 ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_temp.txt > ./BLCA_taxonomy/list_annotated.txt
        source /home/DataAnalysis/miniconda2/bin/activate qiime191
        filter_fasta.py -f DADA2_ASV-seq_fixed.fasta -o ./BLCA_taxonomy/list_unannotated_seqs.fasta -s ./BLCA_taxonomy/list_annotated.txt -n   #extract unannotated sequences
        source /home/DataAnalysis/miniconda2/bin/deactivate qiime191
        python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i ./BLCA_taxonomy/list_unannotated_seqs.fasta \
						           	    --cvrset 0.99 --iset 99 --proc $THREADS \
                                                                    -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
                                                                    -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
                                                                    -o ./BLCA_taxonomy/list_unannotated_seqs_BLCA.txt
        cat ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_temp.txt ./BLCA_taxonomy/list_unannotated_seqs_BLCA.txt > ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment.txt
	rm ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_temp.txt
	rm ./BLCA_taxonomy/list*
done
