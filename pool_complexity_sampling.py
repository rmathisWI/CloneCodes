import random
import csv
from operator import itemgetter
def get_index_number(Seq):
    # """Turn the sequence into a binary string coding for the nucleotides."""
    Dictionary = {'T':'1','C':'2','G':'3','A':'4'}
    return (''.join([Dictionary[Basepair] for Basepair in Seq]))

def compareseq(TestSequence,ReferenceSequence):
    mismatch = 0
    for (n,Base) in enumerate(ReferenceSequence):
        if Base != TestSequence[n]:
            mismatch = mismatch + 1
    return mismatch

f = open('filteredoutput112413TTCACA_ACGGGA.txt')
inputdata = [line.split(',') for line in f]
inputdata.pop(0)

#choose a sequence and find things <m mismatches away
m=8
refseq = inputdata[4828][0]
result1=[]
for bc in inputdata:
    if compareseq(bc[0],refseq)<m:
        result1.append(bc)

for x in range(len(result1)):
    result1[x][1] = int(result1[x][1].strip(' \n'))
    result1[x][0] = [int(g) for g in list(get_index_number(result1[x][0]))]
    #result1[x][0] = [int(g) for g in list(result1[x][0])]
result= []
for barcode in result1:
    temp = []
    for d in barcode[0]:
        temp.append(d)
    temp.append(barcode[1])
    result.append(temp)

result = sorted(result,key=itemgetter(-1),reverse=True)
def WriteToOutput(header,result,outputfilename):
    outputfile = open(outputfilename,'wb')
    resultwriter = csv.writer(outputfile, dialect='excel')
    resultwriter.writerow(header)
    for column in result:
        resultwriter.writerow(column)
    outputfile.close()

WriteToOutput(['Bc','reads'],result,'TTCACA_ACGGGA_filteredoutput.csv')
