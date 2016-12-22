import csv
import random
from numpy import average
from math import log

def ImportDataset(inputfilename):
    file = open(inputfilename)
    csviterator = csv.reader(file, dialect = 'excel')
    result = []
    for row in csviterator:
        result.append(row)
    file.close()
    return result

def ImportFlowCSV(inputfilename,column):
    imported = [x[column] for x in ImportDataset(inputfilename)]
    #log transform
    return [log(float(x),10) for x in imported[1:]]

def SCCFracEp(inputfilename,parE,parM):
    DAPI = ImportFlowCSV(inputfilename,0)
    K818 = ImportFlowCSV(inputfilename,1)
    #Split into DAPI+,-
    D_thresh = 8
    DAPIpos = [K818[n] for n in len(DAPI) if DAPI[n] > D_thresh].sort()
    DAPIneg = [K818[n] for n in len(DAPI) if DAPI[n] > D_thresh].sort()
    DAPIpos_rev = sorted(DAPIpos,reverse=True)
    DAPIneg_rev = sorted(DAPIneg,reverse=True)
    #Find APC thresholds for getting parental E and M sort fractions
    M_gate = DAPIpos[round(len(DAPIpos)*parM)]
    E_gate = DAPIpos_rev[round(len(DAPIpos)*parE)]

    #Use thresholds to find E and M in SCC
    M = 1-(DAPIneg_rev.index(M_gate)+1)/len(DAPIneg)
    E = 1-(DAPIneg.index(E_gate)+1)/len(DAPIneg)

    #Return fraction epithelial
    return E/(E+M)
    
    
    




def MixPopulations(file1,control1,file2,control2,pop_ratio):
    pop1 = [x for x in ImportFlowCSV(file1,1)]
    pop2 = [x for x in ImportFlowCSV(file2,1)]
    control1 = [x for x in ImportFlowCSV(control1,1)]
    control2 = [x for x in ImportFlowCSV(control2,1)]
    control1_avg = average(control1)
    control2_avg = average(control2)

    pop1_norm = [x - control1_avg for x in pop1]
    pop2_norm = [x - control2_avg for x in pop2]

    population_size = min([len(pop1),len(pop2),len(control1),len(control2)])
    
    result = random.sample(pop1_norm,int(population_size * pop_ratio))
    print(len(result))
    [result.append(x) for x in random.sample(pop2_norm,int(population_size * (1 - pop_ratio)))]
    print(len(result))
    return result
    
def ControlFACSdata(file1,control1):
    pop1 = [x for x in ImportFlowCSV(file1,1)]
    control1 = [x for x in ImportFlowCSV(control1,1)]
    control1_avg = average(control1)
    return [x - control1_avg for x in pop1]

def WriteToOutput(result,outputfilename):
    outputfile = open(outputfilename,'w', newline = '')
    resultwriter = csv.writer(outputfile, dialect='excel')
    for row in result:
        resultwriter.writerow([str(row)])
    outputfile.close()
