import sys
seq_fileName = sys.argv[1]
match_score = int(sys.argv[2])
mismatch_score = int(sys.argv[3])
indel_score = int(sys.argv[4])
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
i_max, j_max = len(seq_1) + 10, len(seq_2) + 10;
score_matrix = [[0 for x in range(i_max)] for y in range(j_max)]
score_matrix[0][0] = 0
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
        if (max(cur_score) >= max_score):
            max_score = max(cur_score)
            max_location_i = i
            max_location_j = j
        j += 1
    i += 1

print "Local Alignment Score is: "
#global alignment score is located at score_matrix[len(seq_1)][len(seq_2)]
print max_score
print "Local Alignment ends for seq_1 at index: "
print max_location_i
print "Local Alignment ends for seq_2 at index: "
print max_location_j
