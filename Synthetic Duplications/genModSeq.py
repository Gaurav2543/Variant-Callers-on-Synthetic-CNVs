import os
import random
import numpy as np

def generate_random_dna(length):
    return ''.join([random.choice('ACGT') for _ in range(length)])

# define a function taht converts the sequence to a fasta format and stores it ina file
def seq_to_fasta(header, seq, filename):
    seq = [seq[i:i+80] for i in range(0, len(seq), 80)]
    seq = '\n'.join(seq)
    with open(filename  + '.fasta', 'w') as f:
        f.write(header + seq)

dna = generate_random_dna(1000000)
header = '>Random_DNA_1Mbp\n'
os.makedirs(os.path.dirname('Homozygous/Raw_Seq_Files_CCS/refRand'), exist_ok = True)
seq_to_fasta(header, dna, 'Homozygous/Raw_Seq_Files_CCS/refRand')
os.makedirs(os.path.dirname('Heterozygous/Raw_Seq_Files_CCS/refRand'), exist_ok = True)
seq_to_fasta(header, dna, 'Heterozygous/Raw_Seq_Files_CCS/refRand')

# n is the number of copies of each duplication length
def add_duplications(dna_sequence, copyNumber):
    count = 0
    sum_SV_lengths = 0  
    duplication_info_homo = []
    duplication_info_hetero = []
    duplication_lengths = []
    modified_sequence = dna_sequence  
    
    len1 = random.randint(50, 100)
    len2 = random.randint(300, 500)
    len3 = random.randint(500, 1000)
    len4 = random.randint(1500, 2000)
    len5 = random.randint(4000, 5000)
    len6 = random.randint(8000, 10000)
    len7 = random.randint(25000, 30000)
    
    # copyNumber = 4
    for i in range(copyNumber):
        for j in range(1,8):
            duplication_lengths.append(eval(f'len{j}'))

    num_duplications = 7*copyNumber

    # duplication_positions is a list of 1 Million/(n*7) positions
    duplication_positions = [int(i*len(dna_sequence)/(num_duplications+1)) for i in range(1, num_duplications+1)]

    for duplication_length in duplication_lengths:
        duplicationPosition = duplication_positions[count]
        duplication_end = duplicationPosition + duplication_length
        sequence = modified_sequence[duplicationPosition:duplication_end]
        modified_sequence = modified_sequence[:duplicationPosition] + sequence + modified_sequence[duplicationPosition:]
        sum_SV_lengths += duplication_length
        count += 1
        duplicationPosition += sum_SV_lengths
        duplication_end += sum_SV_lengths
        duplication_info_homo.append((duplicationPosition, duplication_end, duplication_length))
        duplication_info_homo.append((duplicationPosition+len(dna_sequence), duplication_end+len(dna_sequence), duplication_length))
        duplication_info_hetero.append((duplicationPosition, duplication_end, duplication_length))

    return modified_sequence, duplication_info_homo, duplication_info_hetero, sum_SV_lengths
    
modified_dna, duplications_info_homo, duplications_info_hetero, sum_SV_lengths = add_duplications(dna, 5)

# # For Homozygous duplication 
modified_dna_hom = modified_dna + modified_dna
# if directory does not exist, make the direcory using mkdir
seq_to_fasta(header, modified_dna_hom, 'Homozygous/Raw_Seq_Files_CCS/modRand')
print("No. of duplications: ", len(duplications_info_homo))
# for d in duplications_info:
#     print(d)
print("Sum of the lengths of all homozygous structural variations: ", 2*sum_SV_lengths)
print("Length of the modified sequence: ", len(modified_dna_hom))
# Check if the length of the modified sequence is equal to the original sequence minus the sum of the lengths of all structural variations
# print("Homo Duplication Info: ", duplications_info_homo)
print(len(modified_dna_hom)-sum_SV_lengths*2==len(2*dna)) 

# For Heterozygous duplication
modified_dna_het = dna + modified_dna
seq_to_fasta(header, modified_dna_het, 'Heterozygous/Raw_Seq_Files_CCS/modRand')
print("\nNo. of duplications: ", len(duplications_info_hetero))
# for d in duplications_info:
#     print(d)
print("Sum of the lengths of all heterozygous structural variations: ", sum_SV_lengths)
print("Length of the modified sequence: ", len(modified_dna_het))
# Check if the length of the modified sequence is equal to the original sequence minus the sum of the lengths of all structural variations
# print("Hetero Duplication Info: ", duplications_info_hetero)
print(len(modified_dna_het)-sum_SV_lengths==len(2*dna)) 

# For Homozygous duplication 
with open("Homozygous/DUPs.txt", "w") as f:
    f.write("Start\tEnd\tSVLEN\n")
    for i in duplications_info_homo:
        f.write(f"{i[0]}\t{i[1]}\t{i[2]}\n")
    # for i in duplications_info_homo:
    #     f.write(f"{i[0]+len(dna)}\t{i[1]+len(dna)}\t{i[2]}\n")

# For Heterozygous duplication
with open("Heterozygous/DUPs.txt", "w") as f:
    f.write("Start\tEnd\tSVLEN\n")
    for i in duplications_info_hetero:
        f.write(f"{i[0]}\t{i[1]}\t{i[2]}\n")
