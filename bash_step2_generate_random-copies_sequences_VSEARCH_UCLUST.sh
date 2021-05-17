#!/bin/bash

WorkDir=/mnt/home1/qinglong/Test_SimulatedAmplicon_NCBI16SRefSeq_ReadLength_taxonomy
THREADS=20

mkdir $WorkDir/Clustering_Denoising_accuracy

cd $WorkDir

###generate the random abundance for each read
for InputSeqData in $WorkDir/NCBI_16S-rRNA-RefSeq_V1V3.fasta $WorkDir/NCBI_16S-rRNA-RefSeq_V3V5.fasta $WorkDir/NCBI_16S-rRNA-RefSeq_V4.fasta $WorkDir/NCBI_16S-rRNA-RefSeq_V6V9.fasta
do
	InputName=`basename $InputSeqData .fasta`

	while read line1
	read line2
	do
		#generate the random count of each read (maxmium: 50)
		CopyNumber=`echo $(($RANDOM % 50 + 1))`
		#get the name of sequence identifier for each read
		IdentifierName=`echo "$line1" | sed 's/>//g'`
		for count in $(seq 1 $CopyNumber); do echo -e "${line1}-${count}\n$line2"; done > 'Simulated_'$IdentifierName'_sequence.fasta'
		echo -e "$IdentifierName\t$CopyNumber" > 'Simulated_'$IdentifierName'_abundance.txt'
	done < $InputSeqData

	cat Simulated_*_sequence.fasta > $WorkDir/Clustering_Denoising_accuracy/''$InputName'_AbundanceSimulated.fasta' && rm Simulated_*_sequence.fasta
	cat Simulated_*_abundance.txt > $WorkDir/Clustering_Denoising_accuracy/''$InputName'_AbundanceList.txt' && rm Simulated_*_abundance.txt
done


###minicing the output from sequence quality filtering: variable length of amplicon sequences after sliding window-based and maxmium expected error-based filtering
cd $WorkDir/Clustering_Denoising_accuracy

for RegionData in $WorkDir/Clustering_Denoising_accuracy/NCBI_16S-rRNA-RefSeq_V*_AbundanceSimulated.fasta
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
        sed 's/>/>NoTrim-/g' ''$SeqName'_reverse_NoTrim_temp.fasta' > ''$SeqName'_reverse_NoTrim_temp2.fasta'
	/home/qinglong/softwares/bbmap-38.34/reformat.sh in=''$SeqName'_reverse_NoTrim_temp2.fasta' out=''$SeqName'_reverse_NoTrim.fasta' rcomp

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

###since V4 region has 250 nt in length, need to remove trimming length higher than 200 nt for either forward or reverse-complemented reads
for TrimLen in 250 300 350 400 450
do
        rm $WorkDir/Clustering_Denoising_accuracy/NCBI_16S-rRNA-RefSeq_V4_AbundanceSimulated_forward_amplicon/'NCBI_16S-rRNA-RefSeq_V4_forward_'$TrimLen'.fasta'
        rm $WorkDir/Clustering_Denoising_accuracy/NCBI_16S-rRNA-RefSeq_V4_AbundanceSimulated_reverse_amplicon/'NCBI_16S-rRNA-RefSeq_V4_reverse_'$TrimLen'.fasta'
done


###merge the amplicon data with variable read length together for OTU picking
###perform several OTU clustering: VSEARCH-de novo (100% & % 97%), UCLUST-de novo (100% & % 97%), UCLUST-closed-ref (100% & % 97%) and DADA2
for RegionDir in $WorkDir/Clustering_Denoising_accuracy/NCBI_16S-rRNA-RefSeq_*_AbundanceSimulated_*_amplicon
do
	cd $RegionDir
	RegionDirName=`basename $RegionDir`
	cat *.fasta > ''$RegionDirName'_AllLength.fa'

	#"QIIME format": SampleID_SeqID such as PatientID-Day-Collection_SeqID or PatientID.Day.Collection_SeqID
	#"UPARSE/VSEARCH format": SampleID.SeqID such as PatientID_Day_Collection.SeqID
	#In ''$RegionDirName'_AllLength.fa': "100nt-V1V3_NG-041941-100" is fine for QIIME; change to "100nt_V1V3.NG_041941_100" for VSEARCH
	cat ''$RegionDirName'_AllLength.fa' | perl -pe 's/>*(_)+/./g' | perl -pe 's/>*(-)+/_/g' > ''$RegionDirName'_AllLength_VSEARCH.fa'

	###############################################
	#######"VSEARCH_de novo (100% & % 97%)"########
	###############################################
	cd $RegionDir
	for similarity in 0.97 0.98 0.99 1
	do
		mkdir $RegionDir/'VSEARCH2.9_denovo_'$similarity''
		cd $RegionDir/'VSEARCH2.9_denovo_'$similarity''

		#VSEARCH will internally sort the sequences by decreasing the length, cluster at 97% or 100% and relabel with OTU_n
		#generate OTU table and OTU sequences (longest sequence, no size in the identifiers)
		#must use option "--cluster_fast" for input sequences with variable read length from the same region
		#options "--centroids" + "--cluster_fast" will output the longest seed read as OTU sequence
		/home/qinglong/softwares/vsearch-binary-2.9.0/bin/vsearch --threads $THREADS --id $similarity --strand plus --fasta_width 0 \
									  --cluster_fast $RegionDir/''$RegionDirName'_AllLength_VSEARCH.fa' \
									  --relabel OTU_ \
									  --centroids all.otus.fasta \
									  --otutabout all.otutab.txt

	        #assign taxonomy by BLCA with NCBI 16S Microbial DB (downloaded on October 6, 2018)
		mkdir BLCA_taxonomy

		#BLCA current version just use single thread for the script
                #split input into several new files of relatively even size; then run the script on each individual file,
                #your will need to have several bash script for each subseq and later concatenate outputs together
                source /home/DataAnalysis/miniconda2/bin/activate base  #activate conda base environment
                seqkit split --by-part $THREADS all.otus.fasta --out-dir Splits   #split sequences into N parts for downstream parallel processing
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
		rm *.muscle *.fsa *.dblist #remove several temporary files during BLCA taxonomy assignment step if the BLCA haven't remove them

		###unfornatunately, BLCA sometimes fails to annotate some sequences from the input during large run
		cut -d$'\t' -f 1 ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_temp.txt > ./BLCA_taxonomy/list_annotated.txt
		source /home/DataAnalysis/miniconda2/bin/activate qiime191
		filter_fasta.py -f all.otus.fasta -o ./BLCA_taxonomy/list_unannotated_seqs.fasta -s ./BLCA_taxonomy/list_annotated.txt -n   #extract unannotated sequences
		source /home/DataAnalysis/miniconda2/bin/deactivate qiime191
		python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i ./BLCA_taxonomy/list_unannotated_seqs.fasta \
									    --cvrset 0.99 --iset 99 --proc $THREADS \
                                                                            -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
                                                                            -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
                                                                            -o ./BLCA_taxonomy/list_unannotated_seqs_BLCA.txt
		cat ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_temp.txt ./BLCA_taxonomy/list_unannotated_seqs_BLCA.txt > ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment.txt
		rm ./BLCA_taxonomy/list*
		rm ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_temp.txt

		source /home/DataAnalysis/miniconda2/bin/activate qiime191
		#VSEARCH created the OTU table already (all.otutab.txt)
		biom convert -i all.otutab.txt -o otu_table_without_taxa.biom --table-type="OTU table" --to-hdf5

		#make OTU table with BLCA assigned taxonomy
		biom add-metadata -i otu_table_without_taxa.biom \
				  -o otu_table_with_taxa_BLCA.biom \
	     			  --observation-metadata-fp ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment.txt \
				  --observation-header OTUID,taxonomy

		#convert the format of otu table
		biom convert -i otu_table_with_taxa_BLCA.biom -o otu_table_with_taxa_BLCA.txt --to-tsv --header-key taxonomy
		biom convert -i otu_table_without_taxa.biom -o otu_table_without_taxa.txt --to-tsv

		source /home/DataAnalysis/miniconda2/bin/deactivate qiime191
	done

	###############################################
        #######"UCLUST_de novo (100% & % 97%)"#########
        ###############################################
	cd $RegionDir
	source /home/DataAnalysis/miniconda2/bin/activate qiime191  #activate the qiime environment for using UCLUST

        for similarity in 0.97 1
        do
		mkdir $RegionDir/'UCLUST1.2.22_denovo_'$similarity''
                cd $RegionDir/'UCLUST1.2.22_denovo_'$similarity''

		uclust --mergesort $RegionDir/''$RegionDirName'_AllLength.fa' --output $RegionDir/''$RegionDirName'_AllLength_sorted.fa'
		uclust --input $RegionDir/''$RegionDirName'_AllLength_sorted.fa' --uc ''$RegionDirName'_AllLength_sorted.uc' --id $similarity

		#parse the UC file to get the OTU map of each read and each cluster
		python /mnt/home1/qinglong/Test_SimulatedAmplicon_NCBI16SRefSeq_ReadLength_taxonomy/parse_otu_mapping_from_uc.py \
				''$RegionDirName'_AllLength_sorted.uc' ''$RegionDirName'_AllLength_sorted_cluster_map.txt'

		#use longest seed read as representative sequence for OTU
		mkdir rep_set
		pick_rep_set.py -i ''$RegionDirName'_AllLength_sorted_cluster_map.txt' \
				-o ./rep_set/rep_set.fasta \
				-f  $RegionDir/''$RegionDirName'_AllLength_sorted.fa' \
				-m longest
		cat ./rep_set/rep_set.fasta | cut -d ' ' -f 1 > ./rep_set/rep_set_fixed.fasta  #fix the identifier of OTU sequences for "filter_fasta.py"

		rm $RegionDir/''$RegionDirName'_AllLength_sorted.fa'  #remove unnecessary file

                #assign taxonomy by BLCA with NCBI 16S Microbial DB (downloaded on October 6, 2018)
                mkdir BLCA_taxonomy

		#BLCA current version just use single thread for the script
                #split input into several new files of relatively even size; then run the script on each individual file,
                #your will need to have several bash script for each subseq and later concatenate outputs together
                source /home/DataAnalysis/miniconda2/bin/activate base  #activate conda base environment
                seqkit split --by-part $THREADS ./rep_set/rep_set_fixed.fasta --out-dir Splits   #split sequences into N parts for downstream parallel processing
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
                rm *.muscle *.fsa *.dblist #remove several temporary files during BLCA taxonomy assignment step if the BLCA haven't remove them

                ###unfornatunately, BLCA sometimes fails to annotate some sequences from the input during large run
                cut -d$'\t' -f 1 ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_temp.txt > ./BLCA_taxonomy/list_annotated.txt
                source /home/DataAnalysis/miniconda2/bin/activate qiime191
                filter_fasta.py -f ./rep_set/rep_set_fixed.fasta -o ./BLCA_taxonomy/list_unannotated_seqs.fasta -s ./BLCA_taxonomy/list_annotated.txt -n   #extract unannotated sequences
                source /home/DataAnalysis/miniconda2/bin/deactivate qiime191
                python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i ./BLCA_taxonomy/list_unannotated_seqs.fasta \
									    --cvrset 0.99 --iset 99 --proc $THREADS \
                                                                            -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
                                                                            -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
                                                                            -o ./BLCA_taxonomy/list_unannotated_seqs_BLCA.txt
                cat ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_temp.txt ./BLCA_taxonomy/list_unannotated_seqs_BLCA.txt > ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment.txt
                rm ./BLCA_taxonomy/list*
		rm ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_temp.txt

                #create OTU table
            	make_otu_table.py -i ''$RegionDirName'_AllLength_sorted_cluster_map.txt' -o otu_table_without_taxa.biom

                #make OTU table with BLCA assigned taxonomy
                biom add-metadata -i otu_table_without_taxa.biom \
                                  -o otu_table_with_taxa_BLCA.biom \
                                  --observation-metadata-fp ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment.txt \
                                  --observation-header OTUID,taxonomy

                #convert the format of otu table
                biom convert -i otu_table_with_taxa_BLCA.biom -o otu_table_with_taxa_BLCA.txt --to-tsv --header-key taxonomy
		biom convert -i otu_table_without_taxa.biom -o otu_table_without_taxa.txt --to-tsv

	done

	source /home/DataAnalysis/miniconda2/bin/deactivate qiime191

	###############################################
        ######"UCLUST_closed-ref (100% & % 97%)"#######
        ###############################################
        cd $RegionDir
        source /home/DataAnalysis/miniconda2/bin/activate qiime191  #activate the qiime environment for using UCLUST

        for similarity in 0.97 1
        do
                mkdir $RegionDir/'UCLUST1.2.22_ClosedRef_'$similarity''
                cd $RegionDir/'UCLUST1.2.22_ClosedRef_'$similarity''

                uclust --mergesort $RegionDir/''$RegionDirName'_AllLength.fa' --output $RegionDir/''$RegionDirName'_AllLength_sorted.fa'
		#closed-reference picking with SILVA 128 (97_otus_16S.fasta)
                uclust --input $RegionDir/''$RegionDirName'_AllLength_sorted.fa' --uc ''$RegionDirName'_AllLength_sorted.uc' \
			--lib /home/qinglong/16S_RefDB/SILVA128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta \
			--id $similarity \
			--libonly

                #parse the UC file to get the OTU map of each read and each cluster
                python /mnt/home1/qinglong/Test_SimulatedAmplicon_NCBI16SRefSeq_ReadLength_taxonomy/parse_otu_mapping_from_uc_closed-reference.py \
                                ''$RegionDirName'_AllLength_sorted.uc' temp.txt
		cut -d$'\t' -f2- temp.txt > ''$RegionDirName'_AllLength_sorted_cluster_map.txt'  && rm temp.txt

                #use longest seed read as representative sequence for OTU
                mkdir rep_set
                pick_rep_set.py -i ''$RegionDirName'_AllLength_sorted_cluster_map.txt' \
                                -o ./rep_set/rep_set.fasta \
                                -f  $RegionDir/''$RegionDirName'_AllLength_sorted.fa' \
                                -m longest

                rm $RegionDir/''$RegionDirName'_AllLength_sorted.fa'  #remove unnecessary file
		cat ./rep_set/rep_set.fasta | cut -d ' ' -f 1 > ./rep_set/rep_set_fixed.fasta  #fix the identifier of OTU sequences for "filter_fasta.py"

                #assign taxonomy by BLCA with NCBI 16S Microbial DB (downloaded on October 6, 2018)
                mkdir BLCA_taxonomy

		#BLCA current version just use single thread for the script
                #split input into several new files of relatively even size; then run the script on each individual file,
                #your will need to have several bash script for each subseq and later concatenate outputs together
                source /home/DataAnalysis/miniconda2/bin/activate base  #activate conda base environment
                seqkit split --by-part $THREADS ./rep_set/rep_set_fixed.fasta --out-dir Splits   #split sequences into N parts for downstream parallel processing
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
                rm *.muscle *.fsa *.dblist #remove several temporary files during BLCA taxonomy assignment step if the BLCA haven't remove them

                ###unfornatunately, BLCA sometimes fails to annotate some sequences from the input during large run
                cut -d$'\t' -f 1 ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_temp.txt > ./BLCA_taxonomy/list_annotated.txt
                source /home/DataAnalysis/miniconda2/bin/activate qiime191
                filter_fasta.py -f ./rep_set/rep_set_fixed.fasta -o ./BLCA_taxonomy/list_unannotated_seqs.fasta -s ./BLCA_taxonomy/list_annotated.txt -n   #extract unannotated sequences
                source /home/DataAnalysis/miniconda2/bin/deactivate qiime191
                python /home/qinglong/softwares/BLCA-python3/2.blca_main.py -i ./BLCA_taxonomy/list_unannotated_seqs.fasta \
									    --cvrset 0.99 --iset 99 --proc $THREADS \
                                                                            -r /home/qinglong/softwares/BLCA-python3/db/16SMicrobial.ACC.taxonomy \
                                                                            -q /home/qinglong/softwares/BLCA-python3/db/16SMicrobial \
                                                                            -o ./BLCA_taxonomy/list_unannotated_seqs_BLCA.txt
                cat ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_temp.txt ./BLCA_taxonomy/list_unannotated_seqs_BLCA.txt > ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment.txt
                rm ./BLCA_taxonomy/list*
		rm ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_temp.txt

                #create OTU table
                make_otu_table.py -i ''$RegionDirName'_AllLength_sorted_cluster_map.txt' -o otu_table_without_taxa.biom

                #make OTU table with BLCA assigned taxonomy
                biom add-metadata -i otu_table_without_taxa.biom \
                                  -o otu_table_with_taxa_BLCA.biom \
                                  --observation-metadata-fp ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment.txt \
                                  --observation-header OTUID,taxonomy

                #convert the format of otu table
                biom convert -i otu_table_with_taxa_BLCA.biom -o otu_table_with_taxa_BLCA.txt --to-tsv --header-key taxonomy
                biom convert -i otu_table_without_taxa.biom -o otu_table_without_taxa.txt --to-tsv

        done

        source /home/DataAnalysis/miniconda2/bin/deactivate qiime191

done
