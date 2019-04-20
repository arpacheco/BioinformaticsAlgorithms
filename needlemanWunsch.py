#!/usr/bin/env python

""" A simple implementation of the Needleman-Wunsch algorithm for BE 562
course. Accepts two FASTA files as inputs, and outputs sequence alignment plus
scores. """

import sys

base_idx = {'A': 0, 'G': 1, 'C': 2, 'T': 3 }
PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE = 0, 1, 2, 3

def seqalignDP(seq1, seq2, subst_matrix, gap_penalty):
    """
    Return the score of the optimal Needdleman-Wunsch alignment for seq1
    and seq2.
    """
    F = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    TB = [[PTR_NONE for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]

    # initialize dynamic programming table for Needleman-Wunsch alignment
    # (Durbin p.20)
    for i in range(1, len(seq1)+1):
        F[i][0] = 0 - i*gap_penalty
        TB[i][0] = PTR_GAP2  # indicates a gap in seq2
    for j in range(1, len(seq2)+1):
        F[0][j] = 0 - j*gap_penalty
        TB[0][j] = PTR_GAP1  # indicates a gap in seq1

    ## Fill in F and TB Tables
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):

            scores = [0,0,0] #Stores the possible scores from which the maximum will be calculated

            # Match/mismatch
            scores[0] = F[i-1][j-1] + S[base_idx.get(seq1[i-1])][base_idx.get(seq2[j-1])]
            #Gaps
            scores[1] = F[i-1][j] - gap_penalty
            scores[2] = F[i][j-1] - gap_penalty

            F[i][j] = max(scores) #determine the maximum to fill F

            #Fill in the TB matrix
            if max(scores) == scores[0]:
                TB[i][j] = PTR_BASE
            elif max(scores) == scores[1]:
                TB[i][j] = PTR_GAP1
            elif max(scores) == scores[2]:
                TB[i][j] = PTR_GAP2

    return F[len(seq1)][len(seq2)], F, TB

def traceback(seq1, seq2, TB):
    s1 = ""
    s2 = ""

    i = len(seq1)
    j = len(seq2)

    while TB[i][j] != PTR_NONE:
        if TB[i][j] == PTR_BASE:
            s1 = seq1[i-1] + s1
            s2 = seq2[j-1] + s2
            i = i - 1
            j = j - 1
        elif TB[i][j] == PTR_GAP1:
            s1 = '-' + s1
            s2 = seq2[j-1] + s2
            j = j - 1
        elif TB[i][j] == PTR_GAP2:
            s1 = seq1[i-1] + s1
            s2 = '-' + s2
            i = i - 1
        else:
            assert False

    return s1, s2


def readSeq(filename):
    """Reads in a FASTA sequence. Assumes one sequence in the file"""
    seq = []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.rstrip().upper())

    return "".join(seq)

# Substituation matrix and gap_penalty
S = [
    # A   G   C   T
    [ 3, -1, -2, -2],  # A
    [-1,  3, -2, -2],  # G
    [-2, -2,  3, -1],  # C
    [-2, -2, -1,  3]   # T
]
gap_penalty = 4


def main():
    # parse command line
    if len(sys.argv) < 3:
        print("Usage: {0} <FASTA 1> <FASTA 2>".format(sys.argv[0]))
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    seq1 = readSeq(file1)
    seq2 = readSeq(file2)

    score, F, TB = seqalignDP(seq1, seq2, S, gap_penalty)

    s1, s2 = traceback(seq1, seq2, TB)
    print(s1)
    print(s2)

    perfectscore, perfF, perfTB = seqalignDP(seq1, seq1, S, gap_penalty)

    """     To define a distance metric assuming a positive score, the
    perfectscore of a sequence aligned against itself is calculated. The
    actual alignment score was scaled by perfectscore, and this scaled value
    is subtracted from 1. By this algorithm, the distance of a sequence aligned
    against itself is 0.     """

    distance = 1.0-score/float(perfectscore)

    print("Score: {0}".format(score))
    print "perfectscore = "+str(perfectscore)
    print "Distance = "+str(distance)

if __name__ == "__main__":
    main()
