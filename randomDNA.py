import sys
import random
num_sequences = int(sys.argv[1])
size = int(sys.argv[2])
file_1 = open('./randomOut.txt', 'w+')
count_A = 0
count_T = 0
count_C = 0
count_G = 0
file_1.write('\n')
for i in range(0, num_sequences):
    for i in range(0, size):
        num = random.randint(0,99)
        if(num<25):
            file_1.write('A')
            count_A += 1
        if(num>=25 and num<50):
            file_1.write('C')
            count_C += 1
        if(num>=50 and num<75):
            file_1.write('G')
            count_G += 1
        if(num>=75):
            file_1.write('T')
            count_T += 1
    file_1.write('\n')
    file_1.write('\n')
count_A = str(count_A)
count_C = str(count_C)
count_G = str(count_G)
count_T = str(count_T)
file_1.write("Count for each Nucleotide: A[" + count_A + "], C[" + count_C + "], G[" + count_G + "], T[" + count_T + "]" + '\n')
        

