#!/usr/bin/env python

""" This is modified from the bfillings usearch app controller, Qinglong modified for parsing closed-reference clustering by UCLUST

usage: python parse_otu_mapping_from_uc.py X Y

where X is the input .uc file, Y is the output OTU mapping file"""

from sys import argv

def parse_uclust_clusters(clustered_uc_lines):
    """ Returns dict of cluster ID:seq IDs
    clustered_uc_lines: lines from .uc file resulting from closed-reference clustering
    """

    clusters = {}

    seed_hit_ix = 0
    otu_id_ix = 1
    seq_id_ix = 8
    ref_id_ix = 9

    for line in clustered_uc_lines:
        if line.startswith("#") or len(line.strip()) == 0:
            continue
        curr_line = line.strip().split('\t')
        if curr_line[seed_hit_ix] == "L":
            # Need to split on semicolons for sequence IDs to handle case of abundance sorted data
            clusters[curr_line[seq_id_ix]] = [curr_line[seq_id_ix]]
        if curr_line[seed_hit_ix] == "H":
            curr_id = curr_line[seq_id_ix]
            clusters[curr_line[ref_id_ix]].append(curr_id)
    return clusters


input_uc = open(argv[1], "U")
output_mapping = open(argv[2], "w")

clusters = parse_uclust_clusters(input_uc)

for n in clusters:
    curr_seq_ids = "\t".join(clusters[n])
    output_mapping.write("%s\t%s\n" % (n, curr_seq_ids))
