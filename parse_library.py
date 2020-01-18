from __future__ import print_function
from argparse import ArgumentParser
import csv

def parse_primer_sequences(primer_path, name_column, sequence_column):
    primers = []
    with open(primer_path) as primer_file:
        reader = csv.DictReader(primer_file)
        for row in reader:
            primer = {}
            primer['name'] = row[name_column]
            primer['sequence'] = row[sequence_column]
            primers.append(primer)
    return primers

def parse_initiators(initiator_file):
    initiators = []
    with open(initiator_file) as file:
        reader = csv.DictReader(file)
        for row in reader:
            initiators.append([ row['initiator'], row['left sequence'], row['left spacer'], row['right sequence'], row['right spacer'] ])
    return initiators

def parse_sequences(sequence_file):
    with open(sequence_file) as sequences:
        lines = sequences.readlines()
        return [line.strip('\n') for line in lines]

def chop_and_parse_sequences(sequence_file):
    with open(sequence_file) as sequences:
        lines = sequences.readlines()
        sequences = []
        for i in range(0, len(lines)):
            if i % 2 == 0:
                sequences.append(chop_sequence(lines[i], False))
            else:
                sequences.append(chop_sequence(lines[i], True))
        return sequences

def chop_sequence(sequence, probe_first):
    rev_seq = reverseComplement(sequence.strip('\n'))
    chopped_sequence = {}
    chopped_sequence['nt'] = rev_seq[0:20]
    if probe_first:
        chopped_sequence['probe'] = rev_seq[20:45]
        chopped_sequence['initiator'] = rev_seq[47:65]
        chopped_sequence['nb'] = rev_seq[65:83]
    else:
        chopped_sequence['initiator'] = rev_seq[20:38]
        chopped_sequence['probe'] = rev_seq[40:65]
        chopped_sequence['nb'] = rev_seq[65:83]
    return chopped_sequence

def determine_sequence_charachteristics(sequences, pools_and_genes, pools_and_primers, nt_primers, nb_primers, initiators):
    sequence_characteristics = []
    for sequence in sequences:
        initiator = get_initiator(sequence['initiator'], initiators)
        nt_primer = get_primer(sequence['nt'], nt_primers)
        nb_primer = get_primer(sequence['nb'], nb_primers)
        sequence_meta = {}
        sequence_meta['initiator'] = initiator
        sequence_meta['initiator_sequence'] = sequence['initiator']
        sequence_meta['nt_name'] = nt_primer
        sequence_meta['nt_sequence'] = sequence['nt']
        sequence_meta['nb_name'] = nb_primer
        sequence_meta['nb_sequence'] = sequence['nb']
        sequence_meta['probe'] = sequence['probe']
        library_name = get_pool_by_primers(pools_and_primers, nt_primer, nb_primer)
        sequence_meta['gene'] = get_gene_by_pool_and_initiator(pools_and_genes, library_name, initiator)
        sequence_characteristics.append(sequence_meta)
    return sequence_characteristics

def get_initiator(sequence_initiator, initiators):
    for initiator in initiators:
        if sequence_initiator.lower() == initiator[1].lower() or sequence_initiator.lower() == initiator[3].lower():
            return initiator[0]            

def get_primer(sequence_primer, primers):
    for primer in primers:
        if sequence_primer == primer['sequence']:
            return primer['name']

def write_parsed_library(sequences, sequence_meta):
    if len(sequences) != len(sequence_meta):
        raise Exception("Something went wrong")

    with open('parsed_order.csv', 'w') as parsed_file:
        writer = csv.writer(parsed_file)
        writer.writerow(['gene', 'initiator', 'initiator sequence', 'nt name', 'nt sequence', 'nb name', 'nb sequence', 'probe', 'full sequence'])
        for i in range(0, len(sequences)):
            meta = sequence_meta[i]
            sequence = sequences[i]
            writer.writerow([meta['gene'],
                             meta['initiator'], 
                             reverseComplement(meta['initiator_sequence']), 
                             meta['nt_name'],
                             reverseComplement(meta['nt_sequence']), 
                             meta['nb_name'], 
                             reverseComplement(meta['nb_sequence']), 
                             reverseComplement(meta['probe']),
                             sequence])

def get_pool_by_primers(pools_and_primers, nt, nb):
    with open(pools_and_primers) as pools:
        reader = csv.DictReader(pools)
        for row in reader:
            if row['Nb'] == nb and row['Nt'] == nt:
                return row['sublibraries']

def get_gene_by_pool_and_initiator(pool_index, pool, initiator):
    with open(pool_index) as pools:
        reader = csv.DictReader(pools)
        for row in reader:
            if row['name'] == pool:
                return row[initiator.capitalize()]

def reverseComplement(sequence):
    '''
    Get the reverse complement of a sequence. 
    '''
    reverseSequence = reverseString(sequence)
    return getComplement(reverseSequence)

def reverseString(string):
    return string[::-1]

def getComplement(sequence):
    complements = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return ''.join([complements[nuc] for nuc in list(sequence)])

if __name__ == '__main__':
    userInput = ArgumentParser(description="")
    requiredNamed = userInput.add_argument_group('required arguments')
    requiredNamed.add_argument('-s', '--Sequences', action='store', required=True)
    requiredNamed.add_argument('-p', '--Pools', action='store', required=True)
    requiredNamed.add_argument('-g', '--Genes', action='store', required=True)
    requiredNamed.add_argument('-i', '--Initiators', action='store', required=True)
    requiredNamed.add_argument('-nt', '--NtPrimers', action='store', required=True)
    requiredNamed.add_argument('-nb', '--NbPrimers', action='store', required=True)
    args = userInput.parse_args()
    sequence_path = args.Sequences
    pools = args.Pools
    genes = args.Genes
    initiator_path = args.Initiators
    nt_path = args.NtPrimers
    nb_path = args.NbPrimers

    chopped_sequences = chop_and_parse_sequences(sequence_path)
    initiators = parse_initiators(initiator_path)
    nt_primers = parse_primer_sequences(nt_path, 'plate_position', 'sequence')
    nb_primers = parse_primer_sequences(nb_path, 'plate_position', 'sequence')

    sequence_meta = determine_sequence_charachteristics(chopped_sequences, genes, pools, nt_primers, nb_primers, initiators)
    write_parsed_library(parse_sequences(sequence_path), sequence_meta)

## nt:20, initiator:18, probe:25, nb:18
## nt > initiator > probe > nb
## nt > probe > initiator > nb

    