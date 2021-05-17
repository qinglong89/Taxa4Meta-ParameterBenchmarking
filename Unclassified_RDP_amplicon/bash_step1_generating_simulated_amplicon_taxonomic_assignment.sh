#!/bin/bash

WorkDir=/mnt/home1/qinglong/Test_SimulatedAmplicon_NCBI16SRefSeq_ReadLength_taxonomy/Unclassified_amplicon
THREADS=10

cd $WorkDir

#extract uncultured-unclassified bacteria (at least at family level) from RDP 11.5 database
grep '^>' current_Bacteria_unaligned.fa | sed -n '/uncultured/p' | sed -n '/unclassified/p' | cut -d ' ' -f 1 | sed 's/>//g' > family_unclassified_total_identifiers.txt
xargs samtools faidx current_Bacteria_unaligned.fa < family_unclassified_total_identifiers.txt > family_unclassified_total.fasta

#sub-sample 1% of total unclassified bacteria for further testing, but run 10 iterations for such task
source /home/DataAnalysis/miniconda2/bin/activate qiime191
subsample_fasta.py -i family_unclassified_total.fasta -p 0.01 -o family_unclassified.fasta
source /home/DataAnalysis/miniconda2/bin/deactivate qiime191


#manually download the 16S RefSeq data from Bacterial 16S Ribosomal RNA RefSeq Targeted Loci Project (PRJNA33175) and save as "NCBI_16S-rRNA-RefSeq_PRJNA33175_20190709.fasta"
#there are 20,160 sequences in this file, extract sequence identifiers for further mapping
#grep '^>' NCBI_16S-rRNA-RefSeq_PRJNA33175_20190709.fasta | sed 's/>//g' > NCBI_16S-rRNA-RefSeq_PRJNA33175_20190709_SequenceIdentifierList.txt

#>NR_163665.1 Proteus alimentorum strain 08MAS0041 16S ribosomal RNA, partial sequence
#CGAGCGGTAACAGGAGAAAGCTTGCTTTCTTGCTGACGAGCGGCGGACGGGTGAGTAATGTATGGGGATC
#>NR_044847.4 Caedimonas varicaedens strain 221 16S ribosomal RNA, partial sequence
#AGAGTTTGATNNTGGCTCAGAACGAACGCTGGCGGCACGCCTAACACATGCAAGTCGAACGAGGGCATTC
#to format the identiier for OTU picking in QIIME (In this case: change "NR_163665.1 Proteus..." to "100nt-V4_NR-163665" - identifier will be modified following the processes below)
#"QIIME format": SampleID_SeqID such as PatientID-Day-Collection_SeqID or PatientID.Day.Collection_SeqID
#"UPARSE/VSEARCH format": SampleID.SeqID such as PatientID_Day_Collection.SeqID
#cat NCBI_16S-rRNA-RefSeq_PRJNA33175_20190709.fasta | cut -d ' ' -f 1 | perl -pe 's/>*(_)+/-/g' | cut -d '.' -f 1 > NCBI_16S-rRNA-RefSeq.fasta

##########################################################################################################################################################################
###"16S V1-V3 region: forward - AGAGTTTGATCATGGCTCAG, reverse - CCAGCAGCCGCGGTAAT": E. coli 16S V1-V3 amplicon: 490 nt in length
###"16S V3-V5 region: forward - CCTACGGGAGGCAGCAG, reverse - ACTCAAATGAATTGACGG": E. coli 16S V5-V5 amplicon: 551 nt in length
###"16S V4 region: forward - CCAGCAGCCGCGGTAA, reverse - ATTAGATACCCTGGTAGTCC": E. coli 16S V4 amplicon: 253 nt in length
###"16S V6V9 region: forward - AACGCGAAGAACCTTAC, reverse - AAGTCGTAACAAGGTAACCGTA": E. coli 16S V4 amplicon: 507 nt in length
###note that degenerate primers are used for PCR amplification during actual practice, thus need to allow mismatches during pattern search

#use cutadapt to extract the amplicon sequences
cutadapt -a AGAGTTTGATCATGGCTCAG...CCAGCAGCCGCGGTAAT --error-rate 0.2 --cores $THREADS -m 420 -M 510 \
	 --output family-unclassified_V1V3.fasta family-unclassified.fasta
cutadapt -a CCTACGGGAGGCAGCAG...ACTCAAATGAATTGACGG --error-rate 0.2 --cores $THREADS -m 510 -M 561 \
	 --output family-unclassified_V3V5.fasta family-unclassified.fasta
cutadapt -a CCAGCAGCCGCGGTAA...ATTAGATACCCTGGTAGTCC --error-rate 0.2 --cores $THREADS -m 248 -M 258 \
	 --output family-unclassified_V4.fasta family-unclassified.fasta
cutadapt -a AACGCGAAGAACCTTAC...AAGTCGTAACAAGGTAACCGTA --error-rate 0.2 --cores $THREADS -m 477 -M 527 \
	 --output family-unclassified_V6V9.fasta family-unclassified.fasta

#add region to the identifier for further combined OTU picking
sed 's/>/>V1V3_/g' family-unclassified_V1V3.fasta > family-unclassified_V1V3.fa
sed 's/>/>V3V5_/g' family-unclassified_V3V5.fasta > family-unclassified_V3V5.fa
sed 's/>/>V4_/g' family-unclassified_V4.fasta > family-unclassified_V4.fa
sed 's/>/>V6V9_/g' family-unclassified_V6V9.fasta > family-unclassified_V6V9.fa

rm family-unclassified_V*.fasta

##########################################################################################################################################################################
###minicing the output from sequence quality filtering: variable length of amplicon sequences after sliding window-based and maxmium expected error-based filtering
for RegionData in $WorkDir/family-unclassified_V*.fa
do
	SeqName=`basename $RegionData .fa`

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
	rm $WorkDir/family-unclassified_V4_forward_amplicon/'family-unclassified_V4_forward_'$TrimLen'.fasta'
	rm $WorkDir/family-unclassified_V4_reverse_amplicon/'family-unclassified_V4_reverse_'$TrimLen'.fasta'
done

##########################################################################################################################################################################
#########"note: all reverse-complemented reads after length trimming are reverse complemented again to ensure reads with forward orientation for quick alignment"#########
##########################################################################################################################################################################

###perform taxonomic assignment with "BLCA" for each amplicon data
for AmpliconSeqDir in $WorkDir/family-unclassified_*_*_amplicon
do
	cd $AmpliconSeqDir

	for SeqFile in $AmpliconSeqDir/family-unclassified_*_*_*.fasta
	do
	        SeqFileName=`basename $SeqFile .fasta`


	        #"assign taxonomy by BLCA with NCBI 16S Microbial DB (downloaded on July 25, 2019)"
	  	mkdir 'BLCA_taxonomy_'$SeqFileName''

	        python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i $SeqFile \
									    --cvrset 0.95 --iset 97 --proc $THREADS \
		                                                    	    -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
	        	                                           	    -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
	                                                          	    -o ./'BLCA_taxonomy_'$SeqFileName''/''$SeqFileName'_BLCA_out.txt' &  #keep in background

	done

	wait
done

#done for all process, manually combine the BLCA results and RDP taxonomy reults since the sequence identifiers are different across different amplicon data with variable lengths
