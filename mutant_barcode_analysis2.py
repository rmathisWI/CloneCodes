import csv
def WriteToOutput(result,outputfilename):
    outputfile = open(outputfilename,'w',newline='')
    resultwriter = csv.writer(outputfile, dialect='excel')
    for column in result:
        resultwriter.writerow(column)
    outputfile.close()
    

def compareseq(TestSequence,ReferenceSequence):
    return sum([1 for (Tbase,Rbase) in zip(TestSequence,ReferenceSequence) if Tbase != Rbase])
    




sample_number = 13

#Read barcodes 
bcfile = open('barcodes.csv')
csviterator = csv.reader(bcfile, dialect = 'excel')
read_sequences = [[] for samples in range(sample_number)]
for row in csviterator:
    for (i,bc) in enumerate(row):
        read_sequences[i].append(bc)
bcfile.close()

#Read reads for barcodes
bcfile = open('barcodes_reads.csv')
csviterator = csv.reader(bcfile, dialect = 'excel')
read_numbers = [[] for samples in range(sample_number)]
for row in csviterator:
    for (i,bc) in enumerate(row):
        read_numbers[i].append(bc)
bcfile.close()


read_sequences_backup = read_sequences[:]

real_sequences_result = [[] for samples in range(sample_number)]
real_readnumbers_result = [[] for samples in range(sample_number)]

for (sample,sequence_list) in enumerate(read_sequences):
    start_sequences = [seq for seq in sequence_list if len(seq) > 0]
    sequences = start_sequences[:]
    temp_seq_numbers = [int(read_numbers[sample][read_sequences_backup[sample].index(seq)]) for seq in sequences]
    temp_sequences = [seq for seq in sequences]
    for m in range(1,4):
        sequences_result = []
        readnumbers_result = []
        sequences = temp_sequences[:]
        while len(sequences)> 0:
            #Start at the top; that is likely a real sequence
            reference_sequence = sequences[0][:]
            #Find everything similar; this includes the "real" sequence
            mismatches = [read for read in sequences if compareseq(read,reference_sequence) <= m]
            #Now remove all of these identified sequences
            [sequences.remove(mismatch) for mismatch in mismatches]
            #Add the real sequence to the list of real sequences
            sequences_result.append(reference_sequence)
            #Total up the reads from the cluster of similar sequences and add to result list
            readnumbers_result.append(sum([int(temp_seq_numbers[temp_sequences.index(seq)]) for seq in mismatches]))
        temp_seq_numbers = readnumbers_result[:]
        temp_sequences = sequences_result[:]
        print(str(m))
    [real_sequences_result[sample].append(seq) for seq in temp_sequences]
    [real_readnumbers_result[sample].append(reads) for reads in temp_seq_numbers]
    
    print(str(sample))

trans_real_sequences_result = [seq for seq in zip(*real_sequences_result)]
trans_real_readnumbers_result = [seq for seq in zip(*real_readnumbers_result)]

WriteToOutput(trans_real_sequences_result,'nonmut_sequences.csv')
WriteToOutput(trans_real_readnumbers_result,'nonmut_read_numbers.csv')
