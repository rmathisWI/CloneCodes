import csv
def WriteToOutput(result,outputfilename):
    outputfile = open(outputfilename,'w',newline='')
    resultwriter = csv.writer(outputfile, dialect='excel')
    for column in result:
        resultwriter.writerow(column)
    outputfile.close()



notfound = 1E-6
#readsfile = open('reads_for_bc_minus_dropped.csv')
bcfile = open('barcodes.csv')
readsfile = open('barcodes_reads.csv')

csviterator = csv.reader(readsfile, dialect = 'excel')
reads_by_row = []
for row in csviterator:
    reads_by_row.append([number for number in row])
readsfile.close()

all_reads = [[] for samples in range(len(reads_by_row[0]))]
for row in reads_by_row:
    for (i,bc) in enumerate(row):
        all_reads[i].append(bc)


barcodes_by_row = []
csviterator = csv.reader(bcfile, dialect = 'excel')
for row in csviterator:
    barcodes_by_row.append([bc for bc in row])

   
#Find the top barcodes
barcodes = [[] for samples in range(len(barcodes_by_row[0]))]
for row in barcodes_by_row:
    for (i,bc) in enumerate(row):
        barcodes[i].append(bc)
all_barcodes = []
[[all_barcodes.append(bc) for bc in bcs[:]] for bcs in barcodes]
barcodelist500 = list(set(all_barcodes))

#get rid of barcodes that seem to be primer dimers
del barcodelist500[barcodelist500.index('TGGGGGGTGTGACG')]
del barcodelist500[barcodelist500.index('TGTTGTGGGAGTGG')]
del barcodelist500[barcodelist500.index('TTAGGGCCTGTTGA')]

#find reads for barcodes in list
reads = [[] for samples in range(len(barcodelist500))]
for (n,barcode) in enumerate(barcodelist500):
    for i in range(len(all_reads)):
        try:
            reads[n].append(float(all_reads[i][barcodes[i].index(barcode)]))
        except ValueError:
            reads[n].append(1)


#normalize to fraction of total
#calculate totals
total_reads = []
for i in range(len(reads[0])):
    total_reads.append(sum([read[i] for read in reads]))
normread = []
for read in reads:
    normread.append([readnum/total for (readnum,total) in zip(read,total_reads)])


#Now go in and turn the former "1s" into the "not found" fraction.
minimums = [min([read[i] for read in normread]) for i in range(len(normread[0]))]
for (i,reads) in enumerate(normread[:]):
    for (g,read) in enumerate(reads):
        if read == minimums[g]:
            normread[i][g] = notfound

#Get rid of barcodes that are at detection threshold in >=1 tech rep for both cell states at any one day
minimums = [[notfound]*len(normread[0])]
at_minimum = [[mini >= bc for (mini,bc) in zip(minimums[0],barcode)] for barcode in normread]
dropped = 0
for (i,truths) in enumerate(at_minimum):
    for day in [truths[0:4],truths[4:8],truths[8:12]]:
        if day[0:2].count(True) >= 1 and day[2:4].count(True) >= 1:
            normread[i] = [0 for items in normread[0]]
            barcodelist500[i] = '0'
            dropped = dropped + 1
            break
zeros = [0 for thing in range(len(normread[0]))]
while normread.count(zeros) >=1:
    del normread[normread.index(zeros)]
while barcodelist500.count('0') >=1:
    del barcodelist500[barcodelist500.index('0')]


#now normalize to fraction of total again and
#bring to established proportions of total population
population_fractions = [0.5968, 0.5968, 0.4032, 0.4032, 0.5968, 0.5968, 0.4032, 0.4032, 0.5968, 0.5968, 0.4032, 0.4032]
total_reads = [sum([read[i] for read in normread]) for i in range(len(normread[0]))]
norm_filtered_read = []
for read in normread:
    norm_filtered_read.append([readnum/(total/fraction) for (readnum,total,fraction) in zip(read,total_reads,population_fractions)])

WriteToOutput(norm_filtered_read,'barcodes_abundance_prenorm.csv')

#Dealing with contamination
#contamination_fraction = [0.076667704, 0.076667704, 0.07323596, 0.07323596, 0.065371211, 0.065371211, 0.054593875, 0.054593875, 0.081364829, 0.081364829, 0.031791908, 0.031791908]
#contamination_fraction = [0.076667704, 0.07323596, 0.065371211, 0.054593875, 0.081364829, 0.031791908]
#contamination_fraction = [0.0739, 0.0712, 0.065371211, 0.054593875, 0.081364829, 0.031791908]
#contamination_fraction = [0.0771, 0.0637, 0.0512, 0.0605, 0.0479, 0.0553]
contamination_fraction = [0.0648, 0.0488, 0.0441, 0.0263, 0.0353, 0.0260]

#how many contamination reads are there?
totals_reads = [sum([read[i] for read in norm_filtered_read]) for i in range(len(normread[0]))]
#contam_totals = [contam_frac * total_reads for (contam_frac,total_reads) in zip(contamination_fractions,totals_reads)]
contam_reads = []
for read in norm_filtered_read:
    contam_read = []
    contam_read.append(contamination_fraction[0] * (sum(read[2:4])/2))
    contam_read.append(contamination_fraction[0] * (sum(read[2:4])/2))
    contam_read.append(contamination_fraction[1] * (sum(read[0:2])/2))
    contam_read.append(contamination_fraction[1] * (sum(read[0:2])/2))
    contam_read.append(contamination_fraction[2] * (sum(read[6:8])/2))
    contam_read.append(contamination_fraction[2] * (sum(read[6:8])/2))
    contam_read.append(contamination_fraction[3] * (sum(read[4:6])/2))
    contam_read.append(contamination_fraction[3] * (sum(read[4:6])/2))
    contam_read.append(contamination_fraction[4] * (sum(read[10:12])/2))
    contam_read.append(contamination_fraction[4] * (sum(read[10:12])/2))
    contam_read.append(contamination_fraction[5] * (sum(read[8:10])/2))
    contam_read.append(contamination_fraction[5] * (sum(read[8:10])/2))
    contam_reads.append(contam_read)
#subtracting contamination reads
purified_filtered_reads = [[read-contam for (read,contam) in zip(reads,contams)] for (reads,contams) in zip(norm_filtered_read,contam_reads)]



#take negatives to not found threshold
norm_purified_filtered_reads = []

for reads in purified_filtered_reads:
    filtered = reads[:]
    for i in range(len(reads)):
        if reads[i] <= 0:
            filtered[i] = notfound
    norm_purified_filtered_reads.append(filtered)

#normalize to fraction of total
#calculate totals
total_reads = []
for i in range(len(norm_purified_filtered_reads[0])):
    total_reads.append(sum([read[i] for read in norm_purified_filtered_reads]))
final_norm_reads = []
for read in norm_purified_filtered_reads:
    final_norm_reads.append([readnum/(total/fraction) for (readnum,total,fraction) in zip(read,total_reads,population_fractions)])

#bring to established proportions of total population
population_fractions = [0.5968, 0.5968, 0.4032, 0.4032, 0.5968, 0.5968, 0.4032, 0.4032, 0.5968, 0.5968, 0.4032, 0.4032]
totals_reads = [sum([read[i] for read in final_norm_reads]) for i in range(len(normread[0]))]
renorm_read = [[readnum/total for (readnum,total) in zip(read,totals_reads)] for read in final_norm_reads]
final_reads = [[read*fraction for (read,fraction) in zip(reads,population_fractions)] for reads in renorm_read]


#Get rid of barcodes that drop precipitously from d0 to d7 (don't trust them)
dropped2 = 0
sums = [[sum(day) for day in [read[0:4],read[4:8],read[8:12]]] for read in final_reads] 
precipitous_drop = 0.01
for (i,totals) in enumerate(sums[:]):
    try:
        if totals[1]/totals[0] < precipitous_drop:
            dropped = dropped+1
            dropped2 = dropped2+1
            final_reads[i] = [0 for items in final_reads[0]]
            barcodelist500[i] = '0'
    except ZeroDivisionError:
        ignore = 0
        
while final_reads.count(zeros)>= 1:
    del final_reads[final_reads.index(zeros)]
	
while barcodelist500.count('0') >=1:
    del barcodelist500[barcodelist500.index('0')]
print(dropped)



#Now renorm back to population fractions
population_fractions = [0.5968, 0.5968, 0.4032, 0.4032, 0.5968, 0.5968, 0.4032, 0.4032, 0.5968, 0.5968, 0.4032, 0.4032]
totals_reads = [sum([read[i] for read in final_reads]) for i in range(len(normread[0]))]
renorm_reads = [[readnum/total for (readnum,total) in zip(read,totals_reads)] for read in final_reads]
final_renorm_reads = [[read*fraction for (read,fraction) in zip(reads,population_fractions)] for reads in renorm_reads]


#double check totals
totals_reads = [sum([read[i] for read in final_renorm_reads]) for i in range(len(normread[0]))]

#Write reads to output
WriteToOutput(final_renorm_reads,'norm_filtered_read.csv')


#Write barcodes to output
outputfilename = 'norm_filtered_barcodes.csv'
outputfile = open(outputfilename,'w',newline='')
resultwriter = csv.writer(outputfile, dialect='excel')
for column in barcodelist500:
    resultwriter.writerow([str(column)])
outputfile.close()


#looking at stuff lost
minimums = [min([read[i] for read in final_reads]) for i in range(len(normread[0]))]
dropped = []
for i in range(0,len(minimums),2):
    dropped.append([max(read1,read2) for (read1,read2) in zip([read[i] for read in final_reads],[read[i+1] for read in final_reads]) if (read1==minimums[i],read2==minimums[i+1]) in ((True,False),(False,True))])
WriteToOutput(dropped,'dropped.csv')


