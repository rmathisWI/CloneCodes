import datetime

def get_index_number(Seq):
    # """Turn the sequence into a binary string coding for the nucleotides."""
    Dictionary = {'T':'10','C':'01','G':'11','A':'00'}
    return int((''.join([Dictionary[Basepair] for Basepair in Seq])),2)

def deconvolute_index(inputdata,barcodelength):
    #"""Take a number that is an encoded nucleotide sequence (see GetIndexNumber) and translate it back into a sequence"""
    binarysequence = bin(inputdata)[2:]
    binarysequence = (2*barcodelength-len(binarysequence))*'0' + binarysequence #Put 0s back at beginning
    Dictionary = {('1','0'):'T',('0','1'):'C',('1','1'):'G',('0','0'):'A'}
    return ''.join([Dictionary[bincodes] for bincodes in zip(*[iter(binarysequence)]*2)])
      

def analyze_filtered_data(Inputdatafile,OutputDataTitle):
        Possibilities = 1+ int(28*'1',2)
        
        Librarybarcodes = [('AAGACA', 'ACGCCT'), ('AAGACA', 'ACGTCG'), ('ATCACA', 'ACGCGA'), ('CAAACA', 'ACGAGC'), ('CAAACA', 'ACGGCT'), ('CATACA', 'ACGAGG'), ('GTAACA', 'ACGGGA'), ('GTTACA', 'ACGCCT'), ('GTTACA', 'ACGAGG'), ('TAGACA', 'ACGTCC'), ('TTCACA', 'ACGAGC'), ('TTCACA', 'ACGGGA')]

        print(datetime.datetime.now().time())

        for Barcodes in Librarybarcodes:        # Do this for each library barcode set
                f = open(Inputdatafile)
                Possibilities_list = [0 for x in range(Possibilities)]          # Make a list for every possible barcode
                StillMore = 1
                ReadsSoFar = 0
                while StillMore ==1:    
                        Reads = [f.readline()[0:35] for x in range(1000000)]     # Get reads in blocks of 1E6
                        ReadsSoFar = ReadsSoFar + 1000000        
                        if Reads[0] == '':      
                                break
                        Reads = [get_index_number(Read[8:15] + Read[20:27]) for Read in Reads if Read[0:6] == Barcodes[0] and Read[29:35] == Barcodes[1]]
                        #Get the subset of these reads that have the barcode index for this loop, and turn them into the index number
                        for Indexnumber in Reads:               # Add each index number to the list
                                Possibilities_list[Indexnumber] = Possibilities_list[Indexnumber] + 1  
                        if str(ReadsSoFar)[-7:] == '0000000':           
                                print("Reads so far " + str("{:.2e}".format(ReadsSoFar)))
                                print(datetime.datetime.now().time())
                f.close()
                g = open(OutputDataTitle + str(Barcodes[0])+"_" + str(Barcodes[1])+ ".txt",'w')
                g.write(str(Barcodes)+ "\n")
                for (y,bc) in enumerate(Possibilities_list):
                        if bc>1:
                                g.write(str(deconvolute_index(y,14)) + ', ' + str(bc) + "\n")
                        
                print(Barcodes)             
                print(datetime.datetime.now().time())
                

                #Plotted_reads = [x for x in Possibilities_list if x>0]
                #print(str("{:.2e}".format(len(Plotted_reads)))+' Unique Reads' + Barcodes)
                #print(str("{:.2e}".format(sum(Plotted_reads)))+' Total Reads' + Barcodes)

                g.close()

                del Possibilities_list
                
def filter_sequencing_data(Inputdatafile,Outputdatafile,AcceptableScoreThreshold,UnacceptableBases,UnacceptableScoreThreshold):
        print('start')
        print(datetime.datetime.now().time())
        quality_scores_illumina1_7 = "BCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghi"
        AcceptableScores = frozenset(quality_scores_illumina1_7[AcceptableScoreThreshold:])
        UnacceptableScores = frozenset(quality_scores_illumina1_7[:UnacceptableScoreThreshold])

        #First get the data down to smaller reads

        f = open(Inputdatafile,'r')
        g = open(Outputdatafile,'w')
        StillMore = 1
        ReadsSoFar = 0
        FilteredReads = 0
        while StillMore == 1:
                Reads=[]
                ReadsandScores = [f.readline() for x in range(1000000)]
                ReadsSoFar = ReadsSoFar + 1000000
                if ReadsandScores[0] == '':
                        break  #stop when you run out of reads
            #Filter by quality
                for [n,WholeQualityScore] in enumerate(ReadsandScores[3::4]):
                    #if WholeQualityScore == '':
                    #    break
                    QualityScore = WholeQualityScore[0:6] + WholeQualityScore[8:15] + WholeQualityScore[20:27] + WholeQualityScore[29:35]
                    if AcceptableScores.issuperset(QualityScore): 
                        Reads.append(ReadsandScores[4*n+1][0:35])
                    elif sum([QualityScore.count(item) for item in set(QualityScore).difference(AcceptableScores)])<UnacceptableBases and len(set(QualityScore).intersection(UnacceptableScores)) == 0:
                        Reads.append(ReadsandScores[4*n+1][0:35])
                #Filter out reads that don't have the internal barcode or have Ns
                Reads = [Read for Read in Reads if Read[15:20] == 'GAGAG' or Read[15:20] == 'CTCTC' and set(Read).issubset('ATGC') == True]
                for Sequence in Reads:
                        g.write(Sequence + "\n")
                if str(ReadsSoFar)[-7:] == '0000000':
                        print("Lines so far " + str("{:.2e}".format(ReadsSoFar)))
                        print(datetime.datetime.now().time())
                FilteredReads = FilteredReads + len(Reads)
        f.close()
        g.close()
        print("Reads before filtering" + str("{:.2e}".format(ReadsSoFar/4)))
        print("Reads after filtering"  + str("{:.2e}".format(FilteredReads)))

def analyze_sequencing_data(Inputdatafile,Outputdatatitle,AcceptableScoreThreshold,UnacceptableBases,UnacceptableScoreThreshold):
        print("start")
        print(datetime.datetime.now().time())

        Filtered_output_name = "output"+ str(datetime.date.today())+".txt"
        filter_sequencing_data(Inputdatafile,Filtered_output_name,AcceptableScoreThreshold,UnacceptableBases,UnacceptableScoreThreshold)
        print("done filtering")
        print(datetime.datetime.now().time())
        
        analyze_filtered_data(Filtered_output_name,Outputdatatitle)
        print("Done")
        print(datetime.datetime.now().time())
        
