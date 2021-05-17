#!/usr/bin/env python2

# original code source: https://www.biostars.org/p/99886/
# modified by Qinglong Wu for personal use

"""
Convert FASTA to FASTQ file with random quality score from a range (from 30 to 42)
Usage: python fasta2fastq-better.py NAME.fasta NAME.fastq
"""

import sys, os
import random
from Bio import SeqIO

# Get inputs
fa_path = sys.argv[1]
fq_path = sys.argv[2]

# make fastq
with open(fa_path, "r") as fasta, open(fq_path, "w") as fastq:
    for record in SeqIO.parse(fasta, "fasta"):
        # to generate static score (ASCII_BASE=33) within a range (here, from 30 to 42) for each sequence
        #record.letter_annotations["phred_quality"] = [random.randint(30, 42)] * len(record)
        # to generate per-base random score (ASCII_BASE=33) within a range (here, from 30 to 42) for each sequence
        record.letter_annotations["phred_quality"] = [random.randint(30,42) for times in range (len(record))]
        #output old Illumina FASTQ format which encodes PHRED qualities using an ASCII offset of 64
        SeqIO.write(record, fastq, "fastq-illumina")
        #output standard Sanger FASTQ format which encodes PHRED qualities using an ASCII offset of 33
        #SeqIO.write(record, fastq, "fastq-sanger")
