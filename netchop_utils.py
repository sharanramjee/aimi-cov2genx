import logging
import tempfile
import subprocess
from Bio import SeqIO
from collections import defaultdict


def name_peptide(peptide):
    return ''.join([entry[1] for entry in peptide])


def compute_score(peptide):
    l_prob = peptide[0][3]
    r_prob = peptide[-1][3]
    mid_prob = 1
    for entry in peptide[1:-1]:
        mid_prob *= entry[3]
    mid_prob = 1 - mid_prob
    prob = l_prob * mid_prob * r_prob
    return prob


def append_peptide(peptide_dict, genome, end):
    peptide_lens = [8, 9, 10]
    for peptide_len in peptide_lens:
        try:
            peptide = genome[end-peptide_len:end+1]
            if len(peptide) in [9, 10, 11]:
                peptide_dict[name_peptide(peptide[1:])].append(
                    [peptide[-1][4], peptide[1][0], compute_score(peptide), 0.0])
        except IndexError:
            continue
    return peptide_dict


def netchop_parse(seq_ids, netchop_output):
    iter_idx = 0
    line_iterator = iter(netchop_output)
    parsed = list()
    for line in line_iterator:
        if "pos" in line and 'AA' in line and 'score' in line:
            parsed.append([])
            if "----" not in next(line_iterator):
                raise ValueError("Dashes expected")
            line = next(line_iterator)
            while '-------' not in line:
                pos = int(line.split()[0])
                AA = line.split()[1]
                C = line.split()[2]
                score = float(line.split()[3])
                parsed[-1].append([pos, AA, C, score, seq_ids[iter_idx]])
                line = next(line_iterator)
            iter_idx += 1
    return parsed


def form_dict(gene_name):
    seqs = list(SeqIO.parse('genomes/' + gene_name + '_gene_protein.fasta', 'fasta'))
    ids = [seq.id for seq in seqs]
    with open('netchop_out/netchop_output.txt' + gene_name + '_gene') as file:
        netchop_out = file.readlines()
    parsed = netchop_parse(ids, netchop_out)
    peptide_dict = defaultdict(list)
    for genome in parsed:
        for end in range(len(genome)-1, -1, -1):
            if genome[end][2] == 'S':
                peptide_dict = append_peptide(peptide_dict, genome, end)
    return peptide_dict


def netchop_predict(sequences):
    """
    Return netChop predictions for each position in each sequence.
    Parameters
    -----------
    sequences : list of string
        Amino acid sequences to predict cleavage for
    Returns
    -----------
    list of list of float
    The i'th list corresponds to the i'th sequence. Each list gives
    the cleavage probability for each position in the sequence.
    """
    with tempfile.NamedTemporaryFile(suffix=".fsa", mode="w") as input_fd:
        for (i, sequence) in enumerate(sequences):
            input_fd.write("> %d\n" % i)
            input_fd.write(sequence)
            input_fd.write("\n")
        input_fd.flush()
        try:
            output = subprocess.check_output(["netchop", input_fd.name])
        except subprocess.CalledProcessError as e:
            logging.error("Error calling netChop: %s:\n%s" % (e, e.output))
            raise

    parsed = netchop_parse(output)
    assert len(parsed) == len(sequences), \
        "Expected %d results but got %d" % (
            len(sequences), len(parsed))
    assert [len(x) for x in parsed] == [len(x) for x in sequences]
    return parsed

