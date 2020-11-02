#!/usr/bin/env python3

# this is a python script template
# this next line will download the file using curl

gff_file="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
fasta_file="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

import os,gzip,itertools,csv,re

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence



if not os.path.exists(gff_file):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta_file):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")

# Process gff to get gene count, total gene length and coding length
with gzip.open(gff_file, "rt") as fh:
    # now add code to process this
    gff = csv.reader(fh, delimiter="\t")
    gene_count = 0
    gene_total_len = 0
    coding_total_len = 0
    for row in gff:
        if row[0].startswith("#"):
            continue
        if row[2] == 'gene':
            gene_count += 1
            gene_total_len += int(row[4]) - int(row[3])
        if row[2] == 'CDS':
            coding_total_len += int(row[4]) - int(row[3])

# Process fasta to get total genome length
with gzip.open(fasta_file, "rt") as fh:
    pairs = aspairs(fh)
    seqs = dict(pairs)
genome_total_len = len(seqs['Chromosome'])

# Print results
print('Number of genes: {}'.format(gene_count))
print('Total length of the genes: {}'.format(gene_total_len))
print('Total length of genome: {}'.format(genome_total_len))
print('percentage of the genome which is coding: {:.2f}%'.format(100 * (gene_total_len / genome_total_len)))
