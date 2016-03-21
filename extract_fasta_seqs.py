


from sys import stderr, exit, argv
from argparse import ArgumentParser, FileType
import json, itertools
# from Bio import SeqIO
from pyfaidx import Fasta 

in_fasta = 'Homo_sapiens.GRCh38.dna.toplevel.filtered.fa'
hg38 = Fasta(in_fasta)
junction_file = open('junctions.json')
outfile = open('outfile', 'w')


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
def reverse_complement(seq):    
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def extract_junction(junction): 
	s1 = str(hg38[str(j[0])][j[5]:j[6]])
	s2 = str(hg38[str(j[0])][j[7]-1:j[8]])
	if j[1] == '-': 
		s1 = reverse_complement(s1)
		s2 = reverse_complement(s2)
	s = s1 + s2
	return s 


for line in junction_file: 
	j = json.loads(line)		
	s = extract_junction(j)
	outfile.write('>' + j[2] + '_' + j[3] + '_' + j[4] + '\n')
	outfile.write(s + '\n')



print 'done?'














