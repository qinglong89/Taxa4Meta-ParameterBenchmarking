#!/bin/bash
WorkDir=/mnt/home1/qinglong/Test_SimulatedAmplicon_NCBI16SRefSeq_ReadLength_taxonomy/Unclassified_amplicon
THREADS=5

#the only file you need is "current_Bacteria_unaligned.fa" downloaded from RDP release 11.5 database

cd $WorkDir

#extract uncultured-unclassified bacteria (at least at family level) from RDP 11.5 database, this will extract 868,902 sequences that are not classified at least at family level
#grep '^>' current_Bacteria_unaligned.fa | sed -n '/uncultured/p' \
#					| sed -n '/unclassified/p' \
#					| cut -d ' ' -f 1 \
#					| sed 's/>//g' > family_unclassified_total_identifiers.txt

#xargs samtools faidx current_Bacteria_unaligned.fa < family_unclassified_total_identifiers.txt > family_unclassified_total.fasta

#let us perform 10 iterations for the analysis since RDP release 11.5 is too big.
for iteration in 1 2 3 4 5 6 7 8 9 10
do
#        mkdir $WorkDir/'Iteration_'$iteration''
        cd $WorkDir/'Iteration_'$iteration''

	#sub-sample 1% (~ 8,689 sequences) of total unclassified bacteria for further testing, but run 10 iterations for such task
#	source /home/DataAnalysis/miniconda2/bin/activate qiime191
#	subsample_fasta.py -i $WorkDir/family_unclassified_total.fasta -p 0.01 -o family-unclassified.fasta
#	source /home/DataAnalysis/miniconda2/bin/deactivate qiime191

	##########################################################################################################################################################################
	###"16S V1-V3 region: forward - AGAGTTTGATCATGGCTCAG, reverse - CCAGCAGCCGCGGTAAT": E. coli 16S V1-V3 amplicon: 490 nt in length
	###"16S V3-V5 region: forward - CCTACGGGAGGCAGCAG, reverse - ACTCAAATGAATTGACGG": E. coli 16S V5-V5 amplicon: 551 nt in length
	###"16S V4 region: forward - CCAGCAGCCGCGGTAA, reverse - ATTAGATACCCTGGTAGTCC": E. coli 16S V4 amplicon: 253 nt in length
	###"16S V6V9 region: forward - AACGCGAAGAACCTTAC, reverse - AAGTCGTAACAAGGTAACCGTA": E. coli 16S V4 amplicon: 507 nt in length
	###note that degenerate primers are used for PCR amplification during actual practice, thus need to allow mismatches during pattern search

	#use cutadapt to extract the amplicon sequences
#	cutadapt -a AGAGTTTGATCATGGCTCAG...CCAGCAGCCGCGGTAAT --error-rate 0.2 --cores $THREADS -m 420 -M 510 \
#        	 --output family-unclassified_V1V3.fasta family-unclassified.fasta
#	cutadapt -a CCTACGGGAGGCAGCAG...ACTCAAATGAATTGACGG --error-rate 0.2 --cores $THREADS -m 525 -M 561 \
#	         --output family-unclassified_V3V5.fasta family-unclassified.fasta
#	cutadapt -a CCAGCAGCCGCGGTAA...ATTAGATACCCTGGTAGTCC --error-rate 0.2 --cores $THREADS -m 248 -M 258 \
#	         --output family-unclassified_V4.fasta family-unclassified.fasta
#	cutadapt -a AACGCGAAGAACCTTAC...AAGTCGTAACAAGGTAACCGTA --error-rate 0.2 --cores $THREADS -m 477 -M 527 \
#	         --output family-unclassified_V6V9.fasta family-unclassified.fasta

	#add region to the identifier for further combined OTU picking
#	sed 's/>/>V1V3_/g' family-unclassified_V1V3.fasta > family-unclassified_V1V3.fa
#	sed 's/>/>V3V5_/g' family-unclassified_V3V5.fasta > family-unclassified_V3V5.fa
#	sed 's/>/>V4_/g' family-unclassified_V4.fasta > family-unclassified_V4.fa
#	sed 's/>/>V6V9_/g' family-unclassified_V6V9.fasta > family-unclassified_V6V9.fa

#	rm family-unclassified_V*.fasta

	##########################################################################################################################################################################
	cd $WorkDir/'Iteration_'$iteration''

	for AmpliconFile in $PWD/family-unclassified_V*.fa
	do
		AmpliconFileName=`basename $AmpliconFile .fa`

		#first
		mkdir ID90-COV85 && cd ID90-COV85
		cp $AmpliconFile ./
		python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i family-unclassified_V*.fa \
 									    -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
	 								    -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
									    -o ''$AmpliconFileName'_BLCA_out_identity90_cov0.85.txt' \
									    --iset 90 --cvrset 0.85 --proc $THREADS &
		cd ..
		#second
		mkdir ID90-COV90 && cd ID90-COV90
		cp $AmpliconFile ./
        	python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i family-unclassified_V*.fa \
                	        	                                    -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
                        	        	                            -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
                                	        	                    -o ''$AmpliconFileName'_BLCA_out_identity90_cov0.90.txt' \
                                        	        	            --iset 90 --cvrset 0.90 --proc $THREADS &
		cd ..
		#third
		mkdir ID90-COV95 && cd ID90-COV95
		cp $AmpliconFile ./
	        python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i family-unclassified_V*.fa \
	                                                        	    -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
	                               	                                    -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
        	                                                            -o ''$AmpliconFileName'_BLCA_out_identity90_cov0.95.txt' \
                	                               		            --iset 90 --cvrset 0.95 --proc $THREADS &
		cd ..
		#forth
		mkdir ID95-COV95 && cd ID95-COV95
		cp $AmpliconFile ./
	        python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i family-unclassified_V*.fa \
        	                                                	    -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
                	                                               	    -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
                        	                	                    -o ''$AmpliconFileName'_BLCA_out_identity95_cov0.95.txt' \
                                	                	            --iset 95 --cvrset 0.95 --proc $THREADS &
		cd ..
		#fifth
		mkdir ID97-COV95 && cd ID97-COV95
		cp $AmpliconFile ./
	        python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i family-unclassified_V*.fa \
	                                                        	    -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
        	                                     		            -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
                	                                                    -o ''$AmpliconFileName'_BLCA_out_identity97_cov0.95.txt' \
                        	                                            --iset 97 --cvrset 0.95 --proc $THREADS &
               cd ..
               #6th
               mkdir ID97-COV97 && cd ID97-COV97
               cp $AmpliconFile ./
               python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i family-unclassified_V*.fa \
                                                                           -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
                                                                           -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
                                                                           -o ''$AmpliconFileName'_BLCA_out_identity97_cov0.97.txt' \
                                                                           --iset 97 --cvrset 0.97 --proc $THREADS &
	       cd ..
               #7th
               mkdir ID97-COV99 && cd ID97-COV99
               cp $AmpliconFile ./
               python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i family-unclassified_V*.fa \
                                                                           -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
                                                                           -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
                                                                           -o ''$AmpliconFileName'_BLCA_out_identity97_cov0.99.txt' \
                                                                           --iset 97 --cvrset 0.99 --proc $THREADS &
               cd ..
               #8th
               mkdir ID98-COV98 && cd ID98-COV98
               cp $AmpliconFile ./
               python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i family-unclassified_V*.fa \
                                                                           -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
                                                                           -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
                                                                           -o ''$AmpliconFileName'_BLCA_out_identity98_cov0.98.txt' \
                                                                           --iset 98 --cvrset 0.98 --proc $THREADS &
               cd ..
	       #9th
               mkdir ID99-COV99 && cd ID99-COV99
               cp $AmpliconFile ./
               python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i family-unclassified_V*.fa \
                                                                           -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
                                                                           -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
                                                                           -o ''$AmpliconFileName'_BLCA_out_identity99_cov0.99.txt' \
                                                                           --iset 99 --cvrset 0.99 --proc $THREADS &
               #10th
               mkdir ID100-COV100 && cd ID100-COV100
               cp $AmpliconFile ./
               python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i family-unclassified_V*.fa \
                                                                           -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
                                                                           -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
                                                                           -o ''$AmpliconFileName'_BLCA_out_identity100_cov1.00.txt' \
                                                                           --iset 100 --cvrset 1 --proc $THREADS &

               wait

		cd $WorkDir/'Iteration_'$iteration''
		cp ./ID90-COV85/''$AmpliconFileName'_BLCA_out_identity90_cov0.85.txt' ./
		cp ./ID90-COV90/''$AmpliconFileName'_BLCA_out_identity90_cov0.90.txt' ./
		cp ./ID90-COV95/''$AmpliconFileName'_BLCA_out_identity90_cov0.95.txt' ./
		cp ./ID95-COV95/''$AmpliconFileName'_BLCA_out_identity95_cov0.95.txt' ./
		cp ./ID97-COV95/''$AmpliconFileName'_BLCA_out_identity97_cov0.95.txt' ./
		cp ./ID97-COV97/''$AmpliconFileName'_BLCA_out_identity97_cov0.97.txt' ./
                cp ./ID97-COV99/''$AmpliconFileName'_BLCA_out_identity97_cov0.99.txt' ./
                cp ./ID98-COV98/''$AmpliconFileName'_BLCA_out_identity98_cov0.98.txt' ./
                cp ./ID99-COV99/''$AmpliconFileName'_BLCA_out_identity99_cov0.99.txt' ./
                cp ./ID100-COV100/''$AmpliconFileName'_BLCA_out_identity100_cov1.00.txt' ./


		rm -rf ./ID90-COV85/
		rm -rf ./ID90-COV90/
		rm -rf ./ID90-COV95/
		rm -rf ./ID95-COV95/
		rm -rf ./ID97-COV95/
                rm -rf ./ID97-COV97/
                rm -rf ./ID97-COV99/
                rm -rf ./ID98-COV98/
                rm -rf ./ID99-COV99/
		rm -rf ./ID100-COV100/

	done

	#summarize the overclassification rate
	cd $WorkDir/'Iteration_'$iteration''

	for BLCAout in $PWD/*_BLCA_out_identity*_cov*.txt
	do
		BLCAoutName=`basename $BLCAout .txt`
		TotalCount=`cat $BLCAout | wc -l`
		CountUnclassified=`cat $BLCAout | awk -F '\t' '$2~/Unclassified/' | wc -l`
		CountClassified=`cat $BLCAout | sed -n '/Unclassified/!p' | wc -l`
		echo -e "$BLCAoutName\t$CountUnclassified\t$CountClassified\t$TotalCount"
	done > temp
	echo -e "Accession\tUnclassified\tClassified\tTotalCount" | cat - temp > 'Iteration_'$iteration'_overclassification_summary.txt'
	rm temp

done
