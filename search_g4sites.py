import re
import sys

def find_G4_2ndGeneration_pattern(sequence: str):
    pos = 0
    array = []
    g1_pattern = re.compile('GGG+')

    while True:
        m = g1_pattern.search(sequence, pos)
        if m is None:
            break
        else:
            array.append(m.span())
            pos = m.end()
    return array


def find_G4_22pattern_n(sequence: str, n: int):
    pos = 0
    array = []
    g2_pattern = re.compile('GGG+.{1,' + str(n) + '}GGG+')

    while True:
        m = g2_pattern.search(sequence, pos)
        if m is None:
            break
        else:
            array.append(m.span())
            pos = m.end()
    return array


def find_G4_31pattern_n(sequence: str, n: int):
    pos = 0
    array_3 = []
    array_1 = []
    array = []
    g3_pattern = re.compile('GGG+.{1,' + str(n) + '}GGG+.{1,' + str(n) + '}GGG+')
    g1_pattern = re.compile('GGG+')

    pos = 0
    while True:
        m = g1_pattern.search(sequence, pos)
        if m is None:
            break
        else:
            # check the position if it is in G3 motif
            g3_match = g3_pattern.match(sequence, m.start())
            if g3_match is not None:
                array_3.append(g3_match.span())
                array.append(g3_match.span())
                # skip
                pos = g3_match.end()
            else:
                array_1.append(m.span())
                array.append(m.span())
                pos = m.end()
    return (array_3, array_1)


def process_file(fasta_file_path, n: int):
    from Bio import SeqIO
    for record in SeqIO.parse(fasta_file_path, 'fasta'):
        id_part = record.id
        desc_part = record.description
        str_seq = str(record.seq)
        g2_pattern = find_G4_22pattern_n(str_seq, n)
        (g3_pattern, g1_pattern) = find_G4_31pattern_n(str_seq, n)
        print("{}\t{}\t{}\t{}\t{}\t{}".format(id_part, id_part[0:2], desc_part, len(g2_pattern), len(g3_pattern), len(g1_pattern)) )


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type = int, default = 7)
    args, remain = parser.parse_known_args()
    file_list = remain

    #print header 
    print("{}\t{}\t{}\t{}\t{}\t{}".format("id", "type", "description", "G2", "G3", "G1") )
    for filename in file_list:
        process_file(filename, args.n)

