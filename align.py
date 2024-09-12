from os.path import lexists
import requests
import sys

# takes 5 arguments
# (1) File containing the two sequences to be aligned (FASTA format).
# (2) The score for matches;
# (3) The score for mismatches (we will assume that all mismatches are score identically);
# (4) The slippage gap penalty cs;
# (5) The non-slippage gap penalty cn.

# Should print out the optimal alignment score, and the optimal alignment itself

# Inputs:
sequences = list(sys.argv[1:])
match_score = 1
mismatch_score = -1
cs = -1
cn = -2
num_seqs = 0

def get_score(nuc0, nuc1):
    if nuc0 == nuc1:
        return 1
    else:
        return -1

def get_slippage_score(seq, seq_i, DP_val):
    if seq[seq_i-2] == seq[seq_i-1]:
        return cs + DP_val
    else:
        return cn + DP_val

for seq in sequences:
    filename = "sequence" + str(num_seqs) + ".txt"
    # open file if !exists
    if lexists(filename) == False:
        # download file from web
        request = requests.get(seq)
        file = open(filename, 'wb').write(request.content)
        #file.close()
    file = open(filename, 'r')
    # find sequences S and T in the file
    S = file.readlines()[1:2][0].strip('\n')
    file.seek(0)
    T = file.readlines()[4:5][0].strip('\n')
    # run the modified version of the Needleman-Wunsch algorithm
    DP_table = [[0 for i in range(len(S)+1)] for j in range(len(T)+1)]
    DP_table_tracker = [[[-1 for a in range(2)] for b in range(len(S)+1)] for c in range(len(T)+1)]
    for j in range(len(T)+1):
        for i in range(len(S)+1):
            if i > 0 and j > 0:
                # calculate score in relation to upper diagonal
                score_diag = DP_table[j-1][i-1] + get_score(T[j-1], S[i-1])
                # calculate score in relation to above
                score_above = get_slippage_score(T, j, DP_table[j-1][i])
                # calculate score in relation to left
                score_left = get_slippage_score(S, i, DP_table[j][i-1])
                # fill the current DP coordinate
                DP_table[j][i] = max(score_diag, score_above, score_left)
                # track which cell of the DP_table the best score came from
                if DP_table[j][i] == score_diag:
                    DP_table_tracker[j][i][0] = j-1
                    DP_table_tracker[j][i][1] = i-1
                elif DP_table[j][i] == score_above:
                    DP_table_tracker[j][i][0] = j-1
                    DP_table_tracker[j][i][1] = i
                elif DP_table[j][i] == score_left:
                    DP_table_tracker[j][i][0] = j
                    DP_table_tracker[j][i][1] = i-1
            # handle zero cases
            if i == 0:
                DP_table[j][i] = cn*j
            elif j == 0:
                DP_table[j][i] = cn*i
    # give the score of the best alignment of S and T
    final_score = DP_table[len(T)][len(S)]
    print("\nScore of the best alignment of S and T:", final_score)

    # give the best alignment of S and T based on the score and DP_table_tracker
    j = len(T)
    i = len(S)
    # in the worst case the size of the strings will be double the longest one,
    # if every '-' is mapped to some nucleotide
    max_len = max(len(T), len(S))*2
    new_T = ['X' for i in range(max_len)]
    new_S = ['X' for i in range(max_len)]
    S_nucs_processed = 0
    T_nucs_processed = 0
    letter = 0
    while(1):
        # if we arrived at 0,0 in the DP_table, the mapping is complete
        if i == 0 and j == 0:
            break
        # if j will not change then we know to add the next letter from S to new_S
        # and '-' to new_T and vice versa for if i will not change
        next_coordinate = DP_table_tracker[j][i]
        if next_coordinate[0] == j:
            new_T[max_len - letter - 1] = '-'
        else:
            new_T[max_len - letter - 1] = T[len(T) - T_nucs_processed - 1]
            T_nucs_processed = T_nucs_processed + 1
        if next_coordinate[1] == i:
            new_S[max_len - letter - 1] = '-'
        else:
            new_S[max_len - letter - 1] = S[len(S) - S_nucs_processed - 1]
            S_nucs_processed = S_nucs_processed + 1
        letter = letter + 1
        
        # if we reached a border of the DP_table_tracker, then fill in the rest of new_S/T.
        # if we reached the j=0 border then fill the rest of new_S with the rest of S
        # and the rest of new_T with '-'. vice versa if i=0. This completes new_S/T.
        if j == 0:
            for left_in_s in range(i):
                new_S[max_len - letter - left_in_s] = S[len(S) - S_nucs_processed - left_in_s]
                new_T[max_len - letter - left_in_s] = '-'
            break
        if i == 0:
            for left_in_t in range(j):
                new_T[max_len - letter - left_in_t] = T[len(T) - T_nucs_processed - left_in_t]
                new_S[max_len - letter - left_in_t] = '-'
            break

        # keep traversing the DP_table_tracker
        j = next_coordinate[0]
        i = next_coordinate[1]

    # output the results to the user
    S_string = ''.join(new_S).strip('X')
    T_string = ''.join(new_T).strip('X')
    print("The best alignment of S and T is:\nS:", S_string, "\nT:", T_string)
    file.close()
    # repeat above for the following files
    num_seqs = num_seqs + 1
