#!/usr/bin/env python3

import os, gzip, itertools, csv

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

def create_codon_dict():
    codon_dict = dict()
    base = ['A', 'T', 'G', 'C']
    for first in base:
        for second in base:
            for third in base:
                codon=first + second + third
                codon_dict[codon] = 0
    return codon_dict

def genome_description(fasta_file):
    codon_dict = create_codon_dict()
    gene_count = 0
    gene_length = 0
    GC_count = 0
    with gzip.open(fasta_file, "rt") as fh:
        seqs = aspairs(fh)
        for seq in seqs:
            gene_count += 1
            seqname  = seq[0]
            seqstring= seq[1]
            gene_length += len(seqstring)
            for i in range(0, len(seqstring), 3):
                codon = seqstring[i: i + 3]
                codon_dict[codon] += 1
                for nt in codon:
                    if nt == 'G' or nt == 'C':
                        GC_count += 1
    codon_count = sum(codon_dict.values())
    codon_freq = dict()
    for k, v in codon_dict.items():
        codon_freq[k] = v / codon_count

    print('Genome: {}'.format(fasta_file))
    print('Number of genes: {}'.format(gene_count))
    print('Total length of the genes: {}'.format(gene_length))
    print('GC content: {:2f}%'.format(100 * GC_count / gene_length))
    print('Total number of codons : {}'.format(codon_count))
    print()
    return codon_freq

url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))

codon_dict_sp1 = genome_description(file1)
codon_dict_sp2 = genome_description(file2)

with open('codon_frequency.tsv', 'w') as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(('Codon', 'sp1_freq', 'sp2_freq'))
    for key in codon_dict_sp1.keys():
        writer.writerow((key, codon_dict_sp1[key], codon_dict_sp2[key]))
