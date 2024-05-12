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

# n is the number of copies of each deletion length
def add_deletions(dna_sequence, n):
    count = 0
    sum_SV_lengths = 0  
    deletion_info = []
    deletion_lengths = []
    modified_sequence = dna_sequence  
    
    for _ in range(n):
        deletion_lengths.append(random.randint(50, 100))
        deletion_lengths.append(random.randint(300, 500))
        deletion_lengths.append(random.randint(500, 1000))
        deletion_lengths.append(random.randint(1500, 2000))
        deletion_lengths.append(random.randint(4000, 5000))
        deletion_lengths.append(random.randint(8000, 10000))
        deletion_lengths.append(random.randint(25000, 30000))

    num_deletions = 7*n
    # Reduce temp is n is increased
    temp = 750000
    # deletion_positions is a list of 1 Million/(n*7) positions
    deletion_positions = [int(i*temp/(num_deletions+1)) for i in range(1, num_deletions+1)]

    for deletion_length in deletion_lengths:
        deletionPosition = deletion_positions[count]
        if(deletionPosition + deletion_length < len(modified_sequence)):
            deletion_end = deletionPosition + deletion_length
            modified_sequence = modified_sequence[:deletionPosition] + modified_sequence[deletionPosition + deletion_length:]
            sum_SV_lengths += deletion_length
            count += 1
            deletionPosition += sum_SV_lengths
            deletion_end += sum_SV_lengths
            deletion_info.append((deletionPosition, deletion_end, deletion_length))
    
    return modified_sequence, deletion_info, sum_SV_lengths

modified_dna, deletions_info, sum_SV_lengths = add_deletions(dna, 5)

# For Homozygous Deletion 
modified_dna_hom = modified_dna + modified_dna
seq_to_fasta(header, modified_dna_hom, 'Homozygous/Raw_Seq_Files_CCS/modRand')
print("No. of deletions: ", len(deletions_info)+len(deletions_info))
# for d in deletions_info:
#     print(d)
print("Sum of the lengths of all homozygous structural variations: ", 2*sum_SV_lengths)
print("Length of the modified sequence: ", len(modified_dna_hom))
# Check if the length of the modified sequence is equal to the original sequence minus the sum of the lengths of all structural variations
print(sum_SV_lengths*2 + len(modified_dna_hom)==len(2*dna)) 

# For Heterozygous Deletion
modified_dna_het = dna + modified_dna
seq_to_fasta(header, modified_dna_het, 'Heterozygous/Raw_Seq_Files_CCS/modRand')
print("\nNo. of deletions: ", len(deletions_info))
# for d in deletions_info:
#     print(d)
print("Sum of the lengths of all heterozygous structural variations: ", sum_SV_lengths)
print("Length of the modified sequence: ", len(modified_dna_het))
# Check if the length of the modified sequence is equal to the original sequence minus the sum of the lengths of all structural variations
print(sum_SV_lengths + len(modified_dna_het)==len(2*dna)) 

# For Homozygous Deletion 
with open("Homozygous/DELs.txt", "w") as f:
    f.write("Start\tEnd\tSVLEN\n")
    for i in deletions_info:
        f.write(f"{i[0]}\t{i[1]}\t{i[2]}\n")
    for i in deletions_info:
        f.write(f"{i[0]+len(dna)}\t{i[1]+len(dna)}\t{i[2]}\n")

# For Heterozygous Deletion
with open("Heterozygous/DELs.txt", "w") as f:
    f.write("Start\tEnd\tSVLEN\n")
    for i in deletions_info:
        f.write(f"{i[0]}\t{i[1]}\t{i[2]}\n")

