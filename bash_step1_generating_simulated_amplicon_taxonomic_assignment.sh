#!/bin/bash

WorkDir=/mnt/home1/qinglong/Test_SimulatedAmplicon_NCBI16SRefSeq_ReadLength_taxonomy/
THREADS=10

cd $WorkDir

#manually download the 16S RefSeq data from Bacterial 16S Ribosomal RNA RefSeq Targeted Loci Project (PRJNA33175) and save as "NCBI_16S-rRNA-RefSeq_PRJNA33175_20190709.fasta"
#there are 20,160 sequences in this file, extract sequence identifiers for further mapping
grep '^>' NCBI_16S-rRNA-RefSeq_PRJNA33175_20190709.fasta | sed 's/>//g' > NCBI_16S-rRNA-RefSeq_PRJNA33175_20190709_SequenceIdentifierList.txt

#>NR_163665.1 Proteus alimentorum strain 08MAS0041 16S ribosomal RNA, partial sequence
#CGAGCGGTAACAGGAGAAAGCTTGCTTTCTTGCTGACGAGCGGCGGACGGGTGAGTAATGTATGGGGATC
#>NR_044847.4 Caedimonas varicaedens strain 221 16S ribosomal RNA, partial sequence
#AGAGTTTGATNNTGGCTCAGAACGAACGCTGGCGGCACGCCTAACACATGCAAGTCGAACGAGGGCATTC
#to format the identiier for OTU picking in QIIME (In this case: change "NR_163665.1 Proteus..." to "100nt-V4_NR-163665" - identifier will be modified following the processes below)
#"QIIME format": SampleID_SeqID such as PatientID-Day-Collection_SeqID or PatientID.Day.Collection_SeqID
#"UPARSE/VSEARCH format": SampleID.SeqID such as PatientID_Day_Collection.SeqID
cat NCBI_16S-rRNA-RefSeq_PRJNA33175_20190709.fasta | cut -d ' ' -f 1 | perl -pe 's/>*(_)+/-/g' | cut -d '.' -f 1 > NCBI_16S-rRNA-RefSeq.fasta

##########################################################################################################################################################################
###"16S V1-V3 region: forward - AGAGTTTGATCATGGCTCAG, reverse - CCAGCAGCCGCGGTAAT": E. coli 16S V1-V3 amplicon: 490 nt in length
###"16S V3-V5 region: forward - CCTACGGGAGGCAGCAG, reverse - ACTCAAATGAATTGACGG": E. coli 16S V5-V5 amplicon: 551 nt in length
###"16S V4 region: forward - CCAGCAGCCGCGGTAA, reverse - ATTAGATACCCTGGTAGTCC": E. coli 16S V4 amplicon: 253 nt in length
###"16S V6V9 region: forward - AACGCGAAGAACCTTAC, reverse - AAGTCGTAACAAGGTAACCGTA": E. coli 16S V4 amplicon: 507 nt in length
###note that degenerate primers are used for PCR amplification during actual practice, thus need to allow mismatches during pattern search

#use cutadapt to extract the amplicon sequences
cutadapt -a AGAGTTTGATCATGGCTCAG...CCAGCAGCCGCGGTAAT --error-rate 0.2 --cores $THREADS -m 420 -M 510 --output NCBI_16S-rRNA-RefSeq_V1V3_SizeFiltered.fasta NCBI_16S-rRNA-RefSeq.fasta  #keep 90% sequences
cutadapt -a CCTACGGGAGGCAGCAG...ACTCAAATGAATTGACGG --error-rate 0.2 --cores $THREADS -m 525 -M 561 --output NCBI_16S-rRNA-RefSeq_V3V5_SizeFiltered.fasta NCBI_16S-rRNA-RefSeq.fasta  #keep 90% sequences
cutadapt -a CCAGCAGCCGCGGTAA...ATTAGATACCCTGGTAGTCC --error-rate 0.2 --cores $THREADS -m 248 -M 258 --output NCBI_16S-rRNA-RefSeq_V4_SizeFiltered.fasta NCBI_16S-rRNA-RefSeq.fasta  #keep 99% sequences
cutadapt -a AACGCGAAGAACCTTAC...AAGTCGTAACAAGGTAACCGTA --error-rate 0.2 --cores $THREADS -m 477 -M 527 --output NCBI_16S-rRNA-RefSeq_V6V9_SizeFiltered.fasta NCBI_16S-rRNA-RefSeq.fasta #keep 75% sequences

#add region to the identifier for further combined OTU picking
sed 's/>/>V1V3_/g' NCBI_16S-rRNA-RefSeq_V1V3_SizeFiltered.fasta > NCBI_16S-rRNA-RefSeq_V1V3.fasta
sed 's/>/>V3V5_/g' NCBI_16S-rRNA-RefSeq_V3V5_SizeFiltered.fasta > NCBI_16S-rRNA-RefSeq_V3V5.fasta
sed 's/>/>V4_/g' NCBI_16S-rRNA-RefSeq_V4_SizeFiltered.fasta > NCBI_16S-rRNA-RefSeq_V4.fasta
sed 's/>/>V6V9_/g' NCBI_16S-rRNA-RefSeq_V6V9_SizeFiltered.fasta > NCBI_16S-rRNA-RefSeq_V6V9.fasta

rm NCBI_16S-rRNA-RefSeq_*_SizeFiltered.fasta

##########################################################################################################################################################################
###minicing the output from sequence quality filtering: variable length of amplicon sequences after sliding window-based and maxmium expected error-based filtering
for RegionData in $WorkDir/NCBI_16S-rRNA-RefSeq_V*.fasta
do
	SeqName=`basename $RegionData .fasta`

	#Trim the forward sequences, this is typically for "Illumina reads"
	/home/qinglong/softwares/bbmap-38.34/reformat.sh forcetrimright=0 in=$RegionData out=''$SeqName'_forward_NoTrim_temp.fasta'
	sed 's/>/>NoTrim-/g' ''$SeqName'_forward_NoTrim_temp.fasta' > ''$SeqName'_forward_NoTrim.fasta'

	for length1 in 100 150 170 200 250 300 350 400 450
	do
		/home/qinglong/softwares/bbmap-38.34/reformat.sh forcetrimright=$length1 in=''$SeqName'_forward_NoTrim_temp.fasta' out=''$SeqName'_forward_'$length1'_temp.fasta'
		sed 's/>/>'$length1'nt-/g' ''$SeqName'_forward_'$length1'_temp.fasta' > ''$SeqName'_forward_'$length1'.fasta'
	done

	rm *_temp.fasta  #remove temp files first
	mkdir ''$SeqName'_forward_amplicon' && mv *_forward_*.fasta ''$SeqName'_forward_amplicon'

	#Trim the reverse-complemented sequences, this is typically for "454 reads"
	/home/qinglong/softwares/bbmap-38.34/reformat.sh in=$RegionData out=''$SeqName'_reverse_NoTrim_temp.fasta' rcomp    #perform reverse complement
	sed 's/>/>NoTrim-/g' ''$SeqName'_reverse_NoTrim_temp.fasta' > ''$SeqName'_reverse_NoTrim.fasta'

	for length2 in 100 150 170 200 250 300 350 400 450
	do
		/home/qinglong/softwares/bbmap-38.34/reformat.sh forcetrimright=$length2 in=''$SeqName'_reverse_NoTrim_temp.fasta' out=''$SeqName'_reverse_'$length2'_temp.fasta'
		sed 's/>/>'$length2'nt-/g' ''$SeqName'_reverse_'$length2'_temp.fasta' > ''$SeqName'_reverse_'$length2'_temp2.fasta'
		#have to reverse-complemented sequences again as forward sequence for downstream alignment in RDP classifier, BLCA and QIIME ("shorten running time with forward sequence!!!")
		/home/qinglong/softwares/bbmap-38.34/reformat.sh in=''$SeqName'_reverse_'$length2'_temp2.fasta' out=''$SeqName'_reverse_'$length2'.fasta' rcomp  #perform reverse complement
	done

	rm *_temp.fasta *_temp2.fasta  #remove temp files first
	mkdir ''$SeqName'_reverse_amplicon' && mv *_reverse_*.fasta ''$SeqName'_reverse_amplicon'
done

#since V4 region has 250 nt in length, need to remove trimming length higher than 200 nt for either forward or reverse-complemented reads
for TrimLen in 250 300 350 400 450
do
	rm $WorkDir/NCBI_16S-rRNA-RefSeq_V4_forward_amplicon/'NCBI_16S-rRNA-RefSeq_V4_forward_'$TrimLen'.fasta'
	rm $WorkDir/NCBI_16S-rRNA-RefSeq_V4_reverse_amplicon/'NCBI_16S-rRNA-RefSeq_V4_reverse_'$TrimLen'.fasta'
done

##########################################################################################################################################################################
#########"note: all reverse-complemented reads after length trimming are reverse complemented again to ensure reads with forward orientation for quick alignment"#########
##########################################################################################################################################################################

###perform taxonomic assignment with "BLCA" for each amplicon data
for AmpliconSeqDir in $WorkDir/NCBI_16S-rRNA-RefSeq_*_*_amplicon
do
	cd $AmpliconSeqDir

	for SeqFile in $AmpliconSeqDir/NCBI_16S-rRNA-RefSeq_*_*_*.fasta
	do
	        SeqFileName=`basename $SeqFile .fasta`

		#BLCA current version just use single thread for the script
                #split input into several new files of relatively even size; then run the script on each individual file,
                #your will need to have several bash script for each subseq and later concatenate outputs together
                source /home/DataAnalysis/miniconda2/bin/activate base  #activate conda base environment
                seqkit split --by-part 20 $SeqFile --out-dir ''$SeqFileName'_splits'   #split sequences into N parts for downstream parallel processing
                source /home/DataAnalysis/miniconda2/bin/deactivate base

	        #"assign taxonomy by BLCA with NCBI 16S Microbial DB (downloaded on October 6, 2018)"
	  	mkdir 'BLCA_taxonomy_'$SeqFileName''

		for SplitSeq in $AmpliconSeqDir/''$SeqFileName'_splits'/$SeqFileName.part_*.fasta
		do
			SplitSeqName=`basename $SplitSeq .fasta`
	        	python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i $SplitSeq \
		                                                                    -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
	        	                                                            -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
										    --cvrset 0.99 --iset 99 --proc 1 \
	                	                                                    -o ./'BLCA_taxonomy_'$SeqFileName''/''$SplitSeqName'_BLCA_out.txt' &  #keep in background
        	done

		wait  #wait the completion of all the background runs

		cat ./'BLCA_taxonomy_'$SeqFileName''/*_BLCA_out.txt > ./'BLCA_taxonomy_'$SeqFileName''/''$SeqFileName'_BLCA_tax_assignment.txt'
		rm ./'BLCA_taxonomy_'$SeqFileName''/*_BLCA_out.txt  #remove the BLCA results for each sequence split

		rm -rf ''$SeqFileName'_splits'  #remove the sequence splits
		rm *.muscle *.fsa *.dblist #remove several temporary files during BLCA taxonomy assignment step if the BLCA haven't remove them just in case
	done
done

###perform taxonomic assignment with "RDP classifier" for each amplicon data
for AmpliconSeqDir in $WorkDir/NCBI_16S-rRNA-RefSeq_*_*_amplicon
do
        cd $AmpliconSeqDir

        for SeqFile in $AmpliconSeqDir/NCBI_16S-rRNA-RefSeq_*_*_*.fasta
        do
                SeqFileName=`basename $SeqFile .fasta`

		#"assign taxonomy with RDP classifier using SILVA128 DB as the training set"
                source /home/DataAnalysis/miniconda2/bin/activate qiime191   #activate qiime  v1.9.1 version
                #assign_taxonomy.py -m rdp --rdp_max_memory 100000 \
                #                   -i $SeqFile \
                #                   -o 'RDP_taxonomy_'$SeqFileName'' \
                #                   -t /home/qinglong/16S_RefDB/SILVA128_QIIME_release/taxonomy/16S_only/97/consensus_taxonomy_7_levels.txt \
                #                   -r /home/qinglong/16S_RefDB/SILVA128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta
                #above step is very memory intensive and time-consuming process, use alternative command for this step to sppeed up the analysis
                parallel_assign_taxonomy_rdp.py --rdp_max_memory 100000 \
                                                -i $SeqFile \
                                                -o 'RDP_taxonomy_'$SeqFileName'' \
                                                -t /home/qinglong/16S_RefDB/SILVA128_QIIME_release/taxonomy/16S_only/97/consensus_taxonomy_7_levels.txt \
                                                -r /home/qinglong/16S_RefDB/SILVA128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta
	done

	source /home/DataAnalysis/miniconda2/bin/deactivate qiime191 #deactivate qiime environment
done

#done for all process, manually combine the BLCA results and RDP taxonomy reults since the sequence identifiers are different across different amplicon data with variable lengths
