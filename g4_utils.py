import re
from typing import Tuple

def generate_staple_oligomer(
    targetseq_5prime: str, targetseq_3prime: str, linker: str):
    # Generate staple oligomer sequence, 
    # from 5 and 3 prime target site sequence and linker.
    
    from Bio.Seq import Seq
    # in this function, 
    #  variables its name start with seq_ means Seq object.
    seq_5prime = Seq(targetseq_5prime)
    seq_3prime = Seq(targetseq_3prime)

    seq_linker = Seq(linker) 

    ret = seq_3prime.reverse_complement_rna() + seq_linker + seq_5prime.reverse_complement_rna() 
    return str(ret)


def find_ATG(sequence: str):
    pos = 0
    array = []
    atg_pattern = re.compile('ATG')
    while True:
        m = atg_pattern.search(sequence, pos)
        if m is None:
            break
        else:
            array.append(m.span())
            pos = m.end()
    return array


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


def find_G4_22pattern(sequence: str):
    pos = 0
    array = []
    g2_pattern = re.compile('GGG+[ATCGU]{1,7}GGG+')

    while True:
        m = g2_pattern.search(sequence, pos)
        if m is None:
            break
        else:
            array.append(m.span())
            pos = m.end()
    return array


def find_G4_31pattern(sequence: str):
    pos = 0
    array_3 = []
    array_1 = []
    array = []
    g3_pattern = re.compile('GGG+[ATCGU]{1,7}GGG+[ATGCU]{1,7}GGG+')
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


def find_G4_04pattern(sequence: str):
    pass
    pos = 0
    array = []
    g4_pattern = re.compile(
            'GGG+[ATCGU]{1,7}GGG+[ATGCU]{1,7}GGG+[ATGCU]{1,7}GGG+')
    while True:
        m = g4_pattern.search(sequence, pos)
        if m is None:
            break
        else:
            array.append(m.span())
            pos = m.end()
    return array


def get_fragment_downstream(
        seq: str, motif_position: Tuple[int, int],
        length: int = 20, offset: int = 0) -> Tuple[int, int, str]:
    seq_length = len(seq)
    start_pos = motif_position[1] + offset
    end_pos = min(start_pos + length, seq_length)
    return (start_pos, end_pos, seq[start_pos:end_pos])


def get_fragment_upstream(
        seq: str, motif_position: Tuple[int,int],
        length: int = 20, offset: int = 0) -> Tuple[int,int,str]:
    seq_length = len(seq)
    end_pos = motif_position[0] - offset
    start_pos = max(end_pos - length, 0)
    return (start_pos, end_pos, seq[start_pos:end_pos])


def count_GC_content(seq: str) -> int:
    '''
    Count the number of GC base for the given sequence.
    '''
    if seq is None:
        return 0
    seq_upper = seq.upper()
    if not set(seq_upper) <= set('ATCG'):
        if set(seq_upper) <= set('ATGCURYMKSWHBVDN'):
            #XXX For now, ambiguous base has not yet implemented.
            return None
        raise Exception
    return seq_upper.count('C') + seq_upper.count('G')


def calc_RNA_RNA_Stability(
        seq: str, temperature: float = 27.0) -> float:
    '''
    This method computes RNA/RNA stability (deltaG)
    for the given RNA sequence and temperature(Celsius).

    This method assumes that the RNA/RNA are complementary.
    (miss-matches does not exist)
    '''

    RNA_NN_table = {
        "AA": (-6.6, -18.4),
        "AT": (-5.7, -15.5),
        "TA": (-8.1, -22.6),
        "CA": (-10.5, -27.8),
        "CT": (-7.6, -19.2),
        "GA": (-13.3, -35.5),
        "GT": (-10.2, -26.2),
        "CG": (-8.0, -19.4),
        "GC": (-14.2, -34.9),
        "GG": (-12.2, -29.7),

        # Followings are the complementary set of the above NN-pairs.
        "TT": (-6.6, -18.4),
        # "AT" : (-5.7, -15.5),
        # "TA" : (-8.1, -22.6),
        "TG": (-10.5, -27.8),
        "AG": (-7.6, -19.2),
        "TC": (-13.3, -35.5),
        "AC": (-10.2, -26.2),
        # "CG" : (-8.0, -19.4),
        # "GC" : (-14.2, -34.9),
        "CC": (-12.2, -29.7)

    }
    s_init = -10.8
    s_sym = -1.4

    # check the sequence
    seq_upper = seq.upper()
    if not set(seq_upper) <= set('ATCG'):
        if set(seq_upper) <= set('ATGCURYMKSWHBVDN'):
            #XXX For now, ambiguous base has not yet implemented.
            return None
        else:
            raise Exception

    h_delta = float(0.)
    s_delta = float(0.)
    for i in range(len(seq_upper) - 1):
        nn_pair = seq_upper[i:i+2]
        h_delta += RNA_NN_table[nn_pair][0]
        s_delta += RNA_NN_table[nn_pair][1]
    s_delta = s_delta + s_init + s_sym
    g = h_delta - (temperature + 273)*s_delta/1000
    return g


if __name__ == '__main__':
    from Bio import SeqIO
    # fasta_in = 'example_seq/hs_nectin4.fasta'
    fasta_in = 'example_seq/mouse_TRPC6_v1.fasta'
    for record in SeqIO.parse(fasta_in, 'fasta'):
        id_part = record.id
        desc_part = record.description
        seq = record.seq
    str_seq = str(seq)
    g2_pattern = find_G4_22pattern(str_seq)
    (g3_pattern, g1_pattern) = find_G4_31pattern(str_seq)
    #fasta_in = 'example_seq/test_g4.fasta'
    for record in SeqIO.parse(fasta_in, 'fasta'):
        id_part = record.id
        desc_part = record.description
        seq = record.seq
    str_seq = str(seq)
    g4_pattern = find_G4_04pattern(str_seq)
    print(g4_pattern)

    print("G_{} = {} kcal/mol".format(
        25, calc_RNA_RNA_Stability('GCATATGC', 25)))
    print("GC-content of {} {}".format(
        'GCATATGC', count_GC_content('GCATATGC')))
