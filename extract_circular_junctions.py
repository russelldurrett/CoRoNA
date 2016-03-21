

from sys import stderr, exit, argv
from argparse import ArgumentParser, FileType
import json 

gtf_handle = '/Users/red/Desktop/IyerProjects/one+.gtf'
# gtf_handle = '/Users/red/Desktop/IyerProjects/Homo_sapiens.GRCh38.83.gtf'
gtf_file = open(gtf_handle)

trans = {}

# Parse valid exon lines from the GTF file into a dict by transcript_id
for line in gtf_file:
    line = line.strip()
    if not line or line.startswith('#'):
        continue
    if '#' in line:
        line = line.split('#')[0].strip()

    try:
        chrom, source, feature, left, right, score, \
            strand, frame, values = line.split('\t')
    except ValueError:
        continue
    left, right = int(left), int(right)

    if feature != 'exon' or left >= right:
        continue

    values_dict = {}
    for attr in values.split(';')[:-1]:
        attr, _, val = attr.strip().partition(' ')
        values_dict[attr] = val.strip('"')
    # print values_dict

    if 'gene_id' not in values_dict or \
            'transcript_id' not in values_dict:
        continue

    transcript_id = values_dict['transcript_id']
    if transcript_id not in trans:
        trans[transcript_id] = [chrom, strand, [[left, right]]]
    else:
        trans[transcript_id][2].append([left, right])

# Sort exons and merge where separating introns are <=5 bps
for tran, [chrom, strand, exons] in trans.items():
        exons.sort()
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] <= 5:
                tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        trans[tran] = [chrom, strand, tmp_exons]

print 'Transcripts and Exon Coords:'
for t in trans.items(): 
    print t 

# # Calculate and print the unique junctions
exon_exon_junctions = set()
for t, [chrom, strand, texons] in trans.items():
    if strand == '+': 
        for i in range(len(texons)):
            for j in range(i+1, len(texons)):
                exon_exon_junctions.add((chrom, strand, t, 'Exon_' + str(j+1), 'Exon_' + str(i+1), texons[j][0], texons[j][1], texons[i][0], texons[i][1]))
    if strand == '-': 
        for i in range(len(texons)):
            for j in range(i+1, len(texons)):
                exon_exon_junctions.add((chrom, strand, t, 'Exon_' + str(len(texons)-i), 'Exon_' + str(len(texons)-j), texons[i][0], texons[i][1], texons[j][0], texons[j][1]))

exon_exon_junctions = sorted(exon_exon_junctions)

print 'Exon-Exon Junctions:'
for t in exon_exon_junctions: 
    print t



# # Calculate and print the unique junctions
intron_exon_junctions = set()
for t, [chrom, strand, texons] in trans.items():
    if strand == '+': 
        for i in range(len(texons)):
            for j in range(i+1, len(texons)):
                intron_exon_junctions.add((chrom, strand, t, 'Intron_' + str(j), 'Exon_' + str(i+1), texons[j-1][1]+1, texons[j][0]-1, texons[i][0], texons[i][1]))
    if strand == '-': 
        for i in range(len(texons)):
            for j in range(i+1, len(texons)):
                intron_exon_junctions.add((chrom, strand, t, 'Intron_' + str(len(texons)-i-1), 'Exon_' + str(len(texons)-j), texons[i][1]+1, texons[i+1][0]-1, texons[j][0], texons[j][1]))

intron_exon_junctions = sorted(intron_exon_junctions)

print 'Intron-Exon Junctions:'
for t in intron_exon_junctions: 
    print t



# # Calculate and print the unique junctions
exon_intron_junctions = set()
for t, [chrom, strand, texons] in trans.items():
    if strand == '+':
        for i in range(len(texons)):
            for j in range(i+1, len(texons)):
                exon_intron_junctions.add((chrom, strand, t, 'Exon_' + str(j), 'Intron_' + str(i+1), texons[j-1][1]+1, texons[j][0]-1, texons[i][0], texons[i][1]))
    if strand == '-':
        for i in range(len(texons)):
            for j in range(i+1, len(texons)):
                exon_intron_junctions.add((chrom, strand, t, 'Exon_' + str(len(texons)-i-1), 'Intron_' + str(len(texons)-j), texons[i][1]+1, texons[i+1][0]-1, texons[j][0], texons[j][1]))

exon_intron_junctions = sorted(exon_intron_junctions)

print 'Exon-Intron Junctions:'
for t in exon_intron_junctions: 
    print t











