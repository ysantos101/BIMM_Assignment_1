import sys
seq_fileName = sys.argv[1]
match_score = int(sys.argv[2])
mismatch_score = int(sys.argv[3])
indel_score = int(sys.argv[4])
output_align = 0
if(len(sys.argv) > 5):
    if(sys.argv[5] == "-a"):
        output_align = 1
# Add boolean for whether -a is present in args


seq_file = open(seq_fileName, "r")
header = seq_file.readline()
seq_1 = seq_file.readline()
header = seq_file.readline()
header = seq_file.readline()
seq_2 = seq_file.readline()

#print "seq_1 is: "
#print seq_1
#print "\n"
#print "seq_2 is: "
#print seq_2
i_max = len(seq_1) + 10
j_max = len(seq_2) + 10
score_matrix = [[0 for x in range(i_max)] for y in range(j_max)]
score_matrix[0][0] = 0
columns = i_max
rows = j_max
back_matrix = [[0 for foo in range(columns)] for bar in range(rows)]
back_matrix[0][0] = 0
max_score = 0
max_location_i = 0
max_location_j = 0
i = 1
while i <= len(seq_1):
    j = 1
    while j <= len(seq_2):
        cur_score = [0 for m in range(4)]
        cur_score[0] = 0
        if(i-1 >= 0 and j-1 >= 0):
            if(seq_1[i-1] == seq_2[j-1]): # diagonal move, match
                cur_score[1] = score_matrix[i-1][j-1] + match_score
            else:
                cur_score[1] = score_matrix[i-1][j-1] - mismatch_score  #diagonal move, mismatch
        else:
            cur_score[1] = -99999
        if(i-1>= 0):
            cur_score[2] = score_matrix[i-1][j] - indel_score
        else:
            cur_score[2] = -99999
        if(j-1 >= 0):
            cur_score[3] = score_matrix[i][j-1] - indel_score
        else:
            cur_score[3] = -99999
        score_matrix[i][j] = max(cur_score)
        if(max(cur_score) == cur_score[0]):
            back_matrix[i][j] = 4
        if(max(cur_score) == cur_score[1]):
            back_matrix[i][j] = 2
        if(max(cur_score) == cur_score[2]):
            back_matrix[i][j] = 1
        if(max(cur_score) == cur_score[3]):
            back_matrix[i][j] = 3
        if (max(cur_score) >= max_score):
            max_score = max(cur_score)
            max_location_i = i
            max_location_j = j
        j += 1
    i += 1

print "Score for Optimal Local Alignment is: "
#global alignment score is located at score_matrix[len(seq_1)][len(seq_2)]
print max_score
end_location_i = max_location_i
end_location_j = max_location_j

#       ******** Already found score for optimal alignment. 
#                   Already found end_location for optimal local alignment.
#                   Now must trim sequences  and invert matrix, then
#                   repeat above steps to find starting location of optimal
#                   local alignment.


trim_length_seq_1 = len(seq_1) - max_location_i
trim_length_seq_2 = len(seq_2) - max_location_j
seq_1_trim = seq_1[:-trim_length_seq_1]
seq_2_trim = seq_2[:-trim_length_seq_2]

invert_seq_1 = seq_1_trim[::-1]
invert_seq_2 = seq_2_trim[::-1]


i_max = max(len(invert_seq_1), len(invert_seq_2) )
j_max = max(len(invert_seq_2), len(invert_seq_1)) 
reverse_score_matrix = [[0 for x in range(i_max)] for y in range(j_max)]
reverse_score_matrix[0][0] = 0
max_score_reverse = 0
max_location_i = 0
max_location_j = 0
index_i = 0
while index_i < len(invert_seq_1):
    index_j = 0
    while index_j < len(invert_seq_2):
        rev_score = [0 for m in range(4)]
        rev_score[0] = 0
        if(index_i-1 >= 0 and index_j-1 >= 0):
            if(invert_seq_1[index_i-1] == invert_seq_2[index_j-1]): # diagonal move, match
                rev_score[1] = reverse_score_matrix[index_i-1][index_j-1] + match_score
            else:
                rev_score[1] = reverse_score_matrix[index_i-1][index_j-1] - mismatch_score  #diagonal move, mismatch
        else:
            rev_score[1] = -99999
        if(index_i-1 >= 0):
            rev_score[2] = reverse_score_matrix[index_i-1][index_j] - indel_score
        else:
            rev_score[2] = -99999
        if(index_j-1 >= 0):
            rev_score[3] = reverse_score_matrix[index_i][index_j-1] - indel_score
        else:
            rev_score[3] = -99999
        reverse_score_matrix[index_i][index_j] = max(rev_score)
        if (max(rev_score) >= max_score_reverse):
            max_score_reverse = max(rev_score)
            max_location_i = index_i
            max_location_j = index_j
        index_j += 1
    index_i += 1
#The following indices are based off the Inverted, Trimmed strings:

start_location_i = max_location_i
start_location_j = max_location_j

start_location_i = len(seq_1) - trim_length_seq_1 - max_location_i
start_location_j = len(seq_2) - trim_length_seq_2 - max_location_j
#start locations now correspond to forward, trimmed sequences.


#   Backtracking:
rev_alignment_seq_1 = ""
rev_alignment_seq_2 = ""
i = end_location_i -1
j = end_location_j -1

print "Length of the Optimal Local Alignment is: "
align_length = max(end_location_i - start_location_i, end_location_j - start_location_j)
print align_length
if(output_align == 1):
    broken = 0
    while i > start_location_i and i > 0:
        if(broken == 1):
            break
        while j >= start_location_j and j >= 0:
            cur_dir = back_matrix[i][j]
            if(cur_dir == 3):    #move left
                cur_symb_1 = "-"
                cur_symb_2 = seq_2[j]
                j +=  -1
            elif(cur_dir == 1):
                cur_symb_1 = seq_1[i]
                cur_symb_2 = "-"
                i += -1
            elif(cur_dir == 2):
                cur_symb_1 = seq_1[i]
                cur_symb_2 = seq_2[j]
                i += -1
                j += -1
            elif(cur_dir == 4):
                rev_alignment_seq_1 += seq_1[i]
                rev_alignment_seq_2 += seq_2[j]
                broken = 1
                break
                i += -1
                j += -1
            else:
                broken = 1
                break
            rev_alignment_seq_1 += cur_symb_1
            rev_alignment_seq_2 += cur_symb_2
        if(broken == 1):
            break
    alignment_seq_1 =  rev_alignment_seq_1[::-1]
    alignment_seq_2 = rev_alignment_seq_2[::-1]

    print "Alignments are as follows, seq_1 appearing first: "
    print alignment_seq_1
    print alignment_seq_2















































