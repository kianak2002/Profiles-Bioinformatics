import math
import itertools as itl
from itertools import permutations


def get_amino_acids(whole):
    aminos = []
    for i in range(len(whole)):
        if whole[i] not in aminos:
            aminos.append(whole[i])
    if '-' not in aminos:
        aminos.append('-')  # gap is always an answer here:)
    return aminos


def make_profiles(seqs, amino_seqs):
    profiles = {}
    amino_dicts = {}
    for i in range(len(seqs[0])):
        for j in range(len(seqs)):
            if seqs[j][i] not in amino_dicts:
                amino_dicts[seqs[j][i]] = 3  # 1 + pseudo count(2)
            else:
                amino_dicts[seqs[j][i]] += 1

        for m in range(len(amino_seqs)):  # if the frequency of an amino acid is 0
            if amino_seqs[m] not in amino_dicts:
                amino_dicts[amino_seqs[m]] = 2  # pseudo count

        for amino in amino_dicts:
            amino_dicts[amino] = amino_dicts[amino] / (len(seqs) + 2 * len(amino_seqs))
            if amino not in profiles:
                profiles[amino] = []
            profiles[amino].append(amino_dicts[amino])
        amino_dicts = {}

    for amino in profiles:
        sum = 0
        for i in range(len(profiles[amino])):
            sum += profiles[amino][i]
        sum = sum / len(profiles[amino])
        for i in range(len(profiles[amino])):
            profiles[amino][i] = profiles[amino][i] / sum
            profiles[amino][i] = math.log2(profiles[amino][i])

    return profiles


def score(seq, profiles):
    var = 0
    score = 0
    for char in seq:
        score += profiles[char][var]
        var += 1
    return score


def make_all_possible(seq, length):
    seqs_without_gap = []

    res = [seq[i: j] for i in range(len(seq))
           for j in range(i + 1, len(seq) + 1)]

    for sequence in res:
        if len(sequence) <= length:
            
            seqs_without_gap.append(sequence)
    return seqs_without_gap


def add_gap(seqs, length):
    new_seqs = []
    gap = '-'
    for sequence in seqs:
        gap = '-'
        for k in range(length - len(sequence) - 1):
            gap += '-'
        seq_gap = sequence
        if len(sequence) != length:
            seq_gap = sequence + gap
        sth = set(permutations(seq_gap))
        for shit in sth:
            shitty = convert_to_string(shit)
            if check_without_gap(shitty, sequence):
                new_seqs.append(shitty)
    return new_seqs


def convert_to_string(seq):
    new_Seq = ""
    for i in range(len(seq)):
        new_Seq += seq[i]
    return new_Seq


def check_without_gap(seq_gap, seq):
    check = ''
    for char in seq_gap:
        if char != '-':
            check += char
    if check == seq:
        return True
    else:
        return False


def find_max_similarity(seqs, profiles):
    max_score = 0
    best_seq = ""
    for sequence in seqs:
        # print(sequence, ' moshkel mnm')
        score_seq = score(sequence, profiles)
        if score_seq > max_score:
            max_score = score_seq
            best_seq = sequence
    # print(max_score, best_seq)
    return best_seq


if __name__ == '__main__':
    n = int(input())  # number of sequences
    sequences = []
    whole_seqs = ""
    for s in range(n):
        sequences.append(input())  # get sequences
        whole_seqs += sequences[s]
    last_seq = input()
    amino_acids = get_amino_acids(whole_seqs)
    profile = make_profiles(sequences, amino_acids)
    all_seq = make_all_possible(last_seq, len(sequences[0]))
    # print(all_seq)
    all_sequences = add_gap(all_seq, len(sequences[0]))
    # print(all_sequences)
    answer = find_max_similarity(all_sequences, profile)
    print(answer)
