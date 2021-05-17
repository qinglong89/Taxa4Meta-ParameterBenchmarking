#!/bin/bash

WorkDir=/mnt/home1/qinglong/Test_SimulatedAmplicon_NCBI16SRefSeq_ReadLength_taxonomy/Clustering_Denoising_accuracy/NCBI_16S-rRNA-RefSeq_V1V3_AbundanceSimulated_forward_amplicon/

cd $WorkDir

for ClusterDir in $PWD/UCLUST1.2.22_ClosedRef_0.97 $PWD/UCLUST1.2.22_denovo_0.97 $PWD/VSEARCH2.9_denovo_0.99
do
	cd $ClusterDir/BLCA_taxonomy

	sed 's/\t/;/g' otu_seqs_BLCA_tax_assignment.txt > otu_taxonomy.txt	#fix for delimiter for parsing

	#filter each taxa based on confidence score and keep all sorted gene IDs from LCA output for merging purpose
	awk -v OFS="\t" -F ';' '{if ($3 >= 90) print $1,$2; else print $1,"unclassified"}' otu_taxonomy.txt > taxa_domain.txt
	awk -v OFS="\t" -F ';' '{if ($5 >= 90) print $1,$4; else print $1,"unclassified"}' otu_taxonomy.txt > taxa_phylum.txt
	awk -v OFS="\t" -F ';' '{if ($7 >= 90) print $1,$6; else print $1,"unclassified"}' otu_taxonomy.txt > taxa_class.txt
	awk -v OFS="\t" -F ';' '{if ($9 >= 90) print $1,$8; else print $1,"unclassified"}' otu_taxonomy.txt > taxa_order.txt
	awk -v OFS="\t" -F ';' '{if ($11 >= 90) print $1,$10; else print $1,"unclassified"}' otu_taxonomy.txt > taxa_family.txt
	awk -v OFS="\t" -F ';' '{if ($13 >= 90) print $1,$12; else print $1,"unclassified"}' otu_taxonomy.txt > taxa_genus.txt
	awk -v OFS="\t" -F ';' '{if ($15 >= 60) print $1,$14; else print $1,"unclassified"}' otu_taxonomy.txt > taxa_species.txt

	#paste taxa_domain.txt taxa_phylum.txt taxa_class.txt taxa_order.txt taxa_family.txt taxa_genus.txt taxa_species.txt | cut -d$'\t' -f 1 > temp1
	#paste taxa_domain.txt taxa_phylum.txt taxa_class.txt taxa_order.txt taxa_family.txt taxa_genus.txt taxa_species.txt | cut -d$'\t' -f 2,4,6,8,10,12,14 | sed 's/\t/;/g' > temp2
	#paste temp1 temp2 > otu_seqs_BLCA_tax_assignment_parsed.txt

	########################################################################################################################################################################
	#join all ranks together
	#include sort function in join command does not solve this for matching column/field with numeric values, but works for non-nemuric values such as gene1, protein3 etc.
	join -j 1 -a 1 -a 2 -t $'\t' -o auto taxa_domain.txt taxa_phylum.txt | \
        	       join -j 1 -a 1 -a 2 -t $'\t' -o auto - taxa_class.txt | \
        	       join -j 1 -a 1 -a 2 -t $'\t' -o auto - taxa_order.txt | \
        	      join -j 1 -a 1 -a 2 -t $'\t' -o auto - taxa_family.txt | \
        	       join -j 1 -a 1 -a 2 -t $'\t' -o auto - taxa_genus.txt | \
        	     join -j 1 -a 1 -a 2 -t $'\t' -o auto - taxa_species.txt > temp

	cat temp | cut -d$'\t' -f 1 > temp1	#extract OTU ID column

	cat temp | cut -d$'\t' -f 2-8 | sed 's/\t/;/g' > temp2	#extract taxonomy columns and convert to semi-colon delimiter

	paste temp1 temp2 > otu_seqs_BLCA_tax_assignment_parsed.txt
	############################################################################################################################################################################

	rm taxa_*.txt temp* otu_taxonomy.txt

	cd ..

	source /home/DataAnalysis/miniconda2/bin/activate qiime191	#activate qiime1 environment
	#make OTU table with BLCA assigned taxonomy
	biom add-metadata -i otu_table_without_taxa.biom \
        	          -o otu_table_with_taxa_BLCA_parsed.biom \
        	          --observation-metadata-fp ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_parsed.txt \
        	          --observation-header OTUID,taxonomy
	#convert the format of otu table
	biom convert -i otu_table_with_taxa_BLCA_parsed.biom -o otu_table_with_taxa_BLCA_parsed.txt --to-tsv --header-key taxonomy

	#collapse features by their taxonomy at the specified level
	summarize_taxa.py -i otu_table_with_taxa_BLCA_parsed.biom --level 2,3,4,5,6,7 -o ./taxa_summarize
done

