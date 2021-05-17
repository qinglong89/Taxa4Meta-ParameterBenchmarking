#!/bin/bash

WorkDir=/mnt/home1/qinglong/Test_SimulatedAmplicon_NCBI16SRefSeq_ReadLength_taxonomy
THREADS=30

mkdir $WorkDir/Clustering_Denoising_accuracy

cd $WorkDir
###merge the amplicon data with variable read length together for OTU picking
###perform several OTU clustering: VSEARCH-de novo (100% & % 97%), UCLUST-de novo (100% & % 97%), UCLUST-closed-ref (100% & % 97%) and DADA2
for RegionDir in $WorkDir/Clustering_Denoising_accuracy/NCBI_16S-rRNA-RefSeq_*_AbundanceSimulated_*_amplicon
do
        cd $RegionDir
        RegionDirName=`basename $RegionDir`
        ###############################################
        #######"VSEARCH_de novo (100% & % 97%)"########
        ###############################################
        cd $RegionDir
        for similarity in 0.99
        do
                #mkdir $RegionDir/'VSEARCH2.9_denovo_'$similarity''
                cd $RegionDir/'VSEARCH2.9_denovo_'$similarity''

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

                cat ./BLCA_taxonomy/*_BLCA_out.txt > ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment.txt
                rm ./BLCA_taxonomy/*_BLCA_out.txt  #remove the BLCA results for each sequence split

                rm -rf Splits  #remove the sequence splits
        done

        ###############################################
        #######"UCLUST_de novo (100% & % 97%)"#########
        ###############################################
        cd $RegionDir
        source /home/DataAnalysis/miniconda2/bin/activate qiime191  #activate the qiime environment for using UCLUST

        for similarity in 0.97
        do
                #mkdir $RegionDir/'UCLUST1.2.22_denovo_'$similarity''
                cd $RegionDir/'UCLUST1.2.22_denovo_'$similarity''

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

                cat ./BLCA_taxonomy/*_BLCA_out.txt > ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment.txt
                rm ./BLCA_taxonomy/*_BLCA_out.txt  #remove the BLCA results for each sequence split

                rm -rf Splits  #remove the sequence splits
        done

        ###############################################
        ######"UCLUST_closed-ref (100% & % 97%)"#######
        ###############################################
        cd $RegionDir
        source /home/DataAnalysis/miniconda2/bin/activate qiime191  #activate the qiime environment for using UCLUST

        for similarity in 0.97
        do
                #mkdir $RegionDir/'UCLUST1.2.22_ClosedRef_'$similarity''
                cd $RegionDir/'UCLUST1.2.22_ClosedRef_'$similarity''

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

                cat ./BLCA_taxonomy/*_BLCA_out.txt > ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment.txt
                rm ./BLCA_taxonomy/*_BLCA_out.txt  #remove the BLCA results for each sequence split

                rm -rf Splits  #remove the sequence splits
	done

done
