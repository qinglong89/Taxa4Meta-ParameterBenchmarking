#!/usr/bin/env Rscript

#usage: Rscript --vanilla bash_DADA2_script.R /path-to-fastq-files/
#alternative usage: R CMD BATCH bash_DADA2_script.R (but need to specify the path directory in the script)

#define argumens from command line
args <- commandArgs(trailingOnly=TRUE)

#pass the work directory (where fastq files kept) from the command line to the script
#this only take the first argument from the command line
DirPath <- args[1]

#alternative usage (need to deactivate above two commands)
#DirPath <- "/mnt/home1/qinglong/Test_SimulatedAmplicon_NCBI16SRefSeq_ReadLength_taxonomy/Clustering_Denoising_accuracy/NCBI_16S-rRNA-RefSeq_V6V9_AbundanceSimulated_reverse_amplicon"


#this workflow is based on the tutorial of DADA2 version 1.8
library(dada2)
setwd(DirPath)
path <- DirPath
list.files(path)

#get filename
fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))
filenames <- gsub(pattern = "\\.fastq$", "", basename(fnFs))

#set up output for filtered sequence, which will be save in a new dir "dada2filtered"
filtFs <- file.path(path, "dada2Filtered", paste0(filenames, "_filtered.fastq"))
names(filtFs) <- filenames

#filter and trim the sequences (make sure minLen > 5 which is the default k-mers used by DADA2), output is binary file, not traditional fastq file
#note that dada2 does not allow Ns
out <- filterAndTrim(fnFs, filtFs, truncQ=2, maxN=0, maxEE=2, rm.phix=TRUE, multithread=TRUE)

#learn errors
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dev.off()

#dereplicate sequence
derepFs <- derepFastq(filtFs, verbose=TRUE)
names(derepFs) <- filenames

#apply the core sample inference algorithm to the dereplicated data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)	

#construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(dadaFs)
ASVseq <- getSequences(seqtab)
write.table(seqtab, file = "DADA2_ASV-feature-table.txt", sep="\t")
write.table(ASVseq, file = "DADA2_ASV-seq.txt", sep="\t")

quit(save="yes") #this did not change the work directory
