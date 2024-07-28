import re
import sys
import itertools

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


def process_file(fasta_file_path, n: int, filter_func = None):
    from Bio import SeqIO
    ret = []
    for record in SeqIO.parse(fasta_file_path, 'fasta'):
        id_part = record.id
        desc_part = record.description
        str_seq = str(record.seq)
        g2_pattern = find_G4_22pattern_n(str_seq, n)
        (g3_pattern, g1_pattern) = find_G4_31pattern_n(str_seq, n)
        record = {
            'id': id_part, 'description': desc_part,
            "g2": g2_pattern, "g3": g3_pattern, "g1": g1_pattern,
        }
        if filter_func != None:
            if filter_func(record) == True:
                ret.append(record)
            else:
                print("{} omitted".format(record['id']), file = sys.stderr)
        else:
            ret.append(record)
    return ret

def calculate_min_distance(data):
    def calculate_interval(pos1: tuple[int,int], pos2: tuple[int,int]):
        if pos2[0] < pos1[0]:
            pos1, pos2 = pos2, pos1
        return pos2[0] - pos1[1]

    for i in range(len(data)):
        # calculate G2 distance
        record = data[i]
        g2_entries = sorted(record['g2'], key = lambda x: x[0])
        g2_intervals = []
        for g2_first, g2_second in itertools.pairwise(g2_entries):
            g2_intervals.append(g2_second[0] - g2_first[1])
        record["g2_intervals"] = g2_intervals
        record["g2_min_interval"] = min(g2_intervals) if 0 < len(g2_intervals) else None

        g3_entries = sorted(record["g3"], key = lambda x: x[0])
        g1_entries = sorted(record["g1"], key = lambda x: x[0])
        g3_g1_intervals = []
        for g3 in g3_entries:
            for g1 in g1_entries:
                interval = calculate_interval(g3, g1)
                if interval < 0:
                    raise
                g3_g1_intervals.append(interval)
                #print("{} - {} -> {}".format(g3, g1, interval))
        record["g3_g1_min_interval"] = min(g3_g1_intervals) if 0<len(g3_g1_intervals) else None
        
    return data

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type = int, default = 7)
    args, remain = parser.parse_known_args()
    file_list = remain
    print("length of intervals of G-tracts: {}".format(args.n), file = sys.stderr)

    header_entry = [
            "id", "description", "G2", "G3", "G1", "G2_min_interval", "G3_G1_min_interval"
    ]
    print("\t".join(header_entry))

    filter_func = lambda entry: True if entry['id'][0:3] == 'NM_' else False
    for filename in file_list:
        data = process_file(filename, args.n, filter_func)
        calculate_min_distance(data)
        for entry in data:
            record = [
                entry['id'], entry['description'], len(entry['g2']), len(entry['g3']), len(entry['g1']), entry['g2_min_interval'] ,entry['g3_g1_min_interval']
            ]
            entry_str = list(map(str, record))
            s = "\t".join(entry_str)
            print(s)

