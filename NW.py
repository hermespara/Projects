#!/usr/bin/env python

import numpy as np
import argparse



#s1 = "GCATGCT"
#s2 = "GATTACA"
# Arguments for NW :
# change the gap score
# change the match score
# change the mismatch score
parser = argparse.ArgumentParser(description="Needleman Wunsch algorithm")

parser.add_argument("-f1", help = "Join input fasta file 1. If not specified, it will be asked to enter the first sequence.",
                    required = False, type = str)
parser.add_argument("-f2", help = "Join input fasta file 2. If not specified, it will be asked to enter the second sequence.",
                    required = False, type = str)
# Gap penalty options
parser.add_argument("-g", help = "Choose gap penalty (default : -4)",
                    required = False, type = float)

# Gap penalty options
parser.add_argument("-m", help = "Choose match score (default : +2)",
                    required = False, type = float)

parser.add_argument("-mm", help = "Choose mismatch score (default : -3)",
                    required = False, type = float)


args = parser.parse_args()

# Default and optionnal arguments
if args.g:
    gap = args.g
else:
    gap = -4

if args.m:
	match = args.m
else:
	match = 2

if args.mm:
	mismatch = args.mm
else:
	mismatch = -3


# check if two charachters are equal or not
def match_score(alpha, beta):
	if alpha == beta:
		return match
	elif alpha == '-' or beta == '-':
		return gap
	else:
		return mismatch

# prints the alignment
def print_alignment(align1, align2):

    i,j = 0,0

    #calcuate identity, score and aligned sequeces
    symbol = ''
    found = 0
    score = 0
    identity = 0
    for i in range(0,len(align1)):
        # if two AAs are the same, then output "*"
        if align1[i] == align2[i]:
            symbol += '*'
            identity = identity + 1
            score += match_score(align1[i], align2[i])

        # if they are not identical and none of them is gap
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-':
            score += match_score(align1[i], align2[i])
            symbol += ' '
            found = 0

        #if one of them is a gap, output a space
        elif align1[i] == '-' or align2[i] == '-':
            symbol += ' '
            score += gap

    identity = float(identity) / len(align1) * 100

    print('Identity =', "%3.3f" % identity, 'percent')
    print('Score =', score)
    print(align1)
    print(symbol)
    print(align2)

# Implementation of the NW algorithm
def NW(s1,s2):

	N = len(s1)
	M = len(s2)

	Ma = np.zeros(shape=(N+1,M+1)) # Matrix format

	Ma[0,0] = 0
# Matrix initiation
	for j in range(N+1):
		Ma[j,0] = j*gap
	for i in range(M+1):
		Ma[0,i] = i * gap

# Filling the Matrix
	for i in range(1,N+1):
		for j in range(1,M+1):
			left = Ma[i, j-1] + gap
			top = Ma[i-1,j] + gap
			if s1[i-1] == s2[j-1]:
				diag = Ma[i-1, j-1] + match
			else :
				diag = Ma[i-1, j-1] + mismatch

			Ma[i,j] = max(left, top, diag)
# Print the Matrix and the score
	print('Alignment matrix = \n%s\n' % Ma)
	align1, align2 = '', ''
	print('Alignment score : {}'.format(Ma[-1,-1]))

	i,j = N,M
############################
########## Traceback########
############################

	while i > 0 and j > 0:
		score_current = Ma[i][j]
		score_diag = Ma[i-1][j-1]
		score_left = Ma[i][j-1]
		score_up = Ma[i-1][j]

		if score_current == score_diag + match_score(s1[i-1], s2[j-1]):
			a1,a2 = s1[i-1],s2[j-1]
			i,j = i-1,j-1
		elif score_current == score_up + gap:
			a1,a2 = s1[i-1],'-'
			i -= 1
		elif score_current == score_left + gap:
			a1,a2 = '-',s2[j-1]
			j -= 1
		align1 += a1
		align2 += a2


	while i > 0:
		a1,a2 = s1[i-1],'-'
		align1 += a1
		align2 += a2
		i -= 1

	while j > 0:
		a1,a2 = '-',s2[j-1]
		align1 += a1
		align2 += a2
		j -= 1

	print_alignment(align1,align2)


# open fasta file
def fasta_file(fa_file):
	Name = None
	seq = ''
	with open(fa_file, 'r') as file:
		for line in file:
			line = line.strip().upper()
			if line.startswith(">"):
				Name = line[1:]
			else:
				line = line.rstrip()
				seq = seq + line
	return seq

# Main
if args.f1 and args.f2:                     # Manage whether or not there's an input file
    s1= fasta_file(args.f1)
    s2= fasta_file(args.f2)
else:
    s1 = input("Enter sequence 1 :").upper()
    s1=input().upper()
    s2 = input("Enter sequence 2:").upper()

NW(s1,s2)

