#!/usr/bin/env python

""" This is modified from the bfillings usearch app controller

usage: python parse_otu_mapping_from_uc.py X Y

where X is the input .uc file, Y is the output OTU mapping file"""


from sys import argv



def parse_usearch61_clusters(clustered_uc_lines,
                             otu_prefix='denovo',
                             ref_clustered=False):
    """ Returns dict of cluster ID:seq IDs
    clustered_uc_lines: lines from .uc file resulting from de novo clustering
    otu_prefix: string added to beginning of OTU ID.
    ref_clustered: If True, will attempt to create dict keys for clusters as
     they are read from the .uc file, rather than from seed lines.
    """

    clusters = {}
    failures = []

    seed_hit_ix = 0
    otu_id_ix = 1
    seq_id_ix = 8
    ref_id_ix = 9

    for line in clustered_uc_lines:
        if line.startswith("#") or len(line.strip()) == 0:
            continue
        curr_line = line.strip().split('\t')
        if curr_line[seed_hit_ix] == "S":
            # Need to split on semicolons for sequence IDs to handle case of
            # abundance sorted data
            clusters[otu_prefix + curr_line[otu_id_ix]] =\
                [curr_line[seq_id_ix].split(';')[0].split()[0]]
        if curr_line[seed_hit_ix] == "H":
            curr_id = curr_line[seq_id_ix].split(';')[0].split()[0]
            if ref_clustered:
                try:
                    clusters[otu_prefix + curr_line[ref_id_ix]].append(curr_id)
                except KeyError:
                    clusters[otu_prefix + curr_line[ref_id_ix]] = [curr_id]
            else:
                clusters[otu_prefix +
                         curr_line[otu_id_ix]].append(curr_id)
        if curr_line[seed_hit_ix] == "N":
            failures.append(curr_line[seq_id_ix].split(';')[0])

    return clusters, failures


input_uc = open(argv[1], "U")
output_mapping = open(argv[2], "w")

clusters,failures = parse_usearch61_clusters(input_uc)

if(len(failures)>0):
    raise(ValueError,"Failures present\n%s" % failures)

for n in clusters:
    curr_seq_ids = "\t".join(clusters[n])
    output_mapping.write("%s\t%s\n" % (n, curr_seq_ids))
