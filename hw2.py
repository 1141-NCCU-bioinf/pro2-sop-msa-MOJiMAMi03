import numpy as np
import pandas as pd

def calculate_SoP(input_path, score_path, gopen, gextend):
    sequences = read_fasta(input_path)
    score_matrix = read_scoring_matrix(score_path)
    
    total_score = 0
    n_sequences = len(sequences)
    
    for i in range(n_sequences):
        for j in range(i + 1, n_sequences):
            pair_score = calculate_pairwise_score(
                sequences[i], 
                sequences[j], 
                score_matrix, 
                gopen, 
                gextend
            )
            total_score += pair_score
    
    return total_score


def read_fasta(filepath):
    sequences = []
    current_seq = ""
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ""
            else:
                current_seq += line
        
        if current_seq:
            sequences.append(current_seq)
    
    return sequences


def read_scoring_matrix(filepath):
    score_matrix = {}
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    lines = [line.strip() for line in lines if line.strip() and not line.strip().startswith('#')]
    
    header = lines[0].split()
    
    for line in lines[1:]:
        parts = line.split()
        row_aa = parts[0].upper()
        scores = parts[1:]
        
        for j, score_str in enumerate(scores):
            if j < len(header):
                col_aa = header[j].upper()
                try:
                    score = int(score_str)
                    score_matrix[(row_aa, col_aa)] = score
                    score_matrix[(col_aa, row_aa)] = score
                except ValueError:
                    continue
    
    return score_matrix


def calculate_pairwise_score(seq1, seq2, score_matrix, gopen, gextend):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have the same length in MSA")
    
    total_score = 0
    in_gap = False
    
    for k in range(len(seq1)):
        a = seq1[k].upper()
        b = seq2[k].upper()
        
        if a != '-' and b != '-':
            if (a, b) in score_matrix:
                total_score += score_matrix[(a, b)]
            in_gap = False
        else:
            if in_gap:
                total_score += gextend
            else:
                total_score += gopen
                in_gap = True
    
    return total_score

if __name__ in "__main__":
    print(calculate_SoP("examples/test1.fasta", "examples/pam250.txt", -10, -2))