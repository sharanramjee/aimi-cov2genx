import Levenshtein, sys
from calculations.python.paths import *
from calculations.python.util import *

with open(data_path(REFERENCE_GENOME), 'r') as file:
	reference = ''.join(file.read().split('\n')[1:])

with open(data_path(REFERENCE_GENES), 'r') as file:
	genes = [g.split() for g in file.read().split('\n')]
	genes = {g[0]: translate(transcribe(reference[int(g[1]):int(g[2])])) for g in genes}

genes = {k: v for k, v in genes.items() if k in ["S_gene", "E_gene", "M_gene", "N_gene"]}

already_seen = {}

columns = next(sys.stdin)[:-1]
print(columns, end=',edit_distance\n')
columns = columns.split(',')

peptide_column = columns.index("peptide")
gene_column = columns.index("gene")

already_seen = {}
for line in sys.stdin:
	line = line[:-1]
	line_list = line.split(',')
	peptide = line_list[peptide_column]
	gene = line_list[gene_column]
	if (peptide, gene) in already_seen:
		other_d = already_seen[(peptide, gene)]
	else:
		lev = lambda x: Levenshtein.distance(peptide, genes[gene][x: x + len(peptide)])
		other_d = min((pos for pos in range(len(genes[gene]) - len(peptide))), key=lev)
		other_d = lev(other_d)
		already_seen[(peptide, gene)] = other_d
	print(line, end=',')
	print(other_d)
