import random
import csv
from math import pow
from math import log
from itertools import chain

def ImportCloneDataset(inputfilename):
    file = open(inputfilename)
    csviterator = csv.reader(file, dialect = 'excel')
    result = []
    for row in csviterator:
        result.append([float(row[0]),float(row[1]),float(row[2])])
    file.close()
    return result

def WriteToOutput(header,result,outputfilename):
    outputfile = open(outputfilename,'wb')
    resultwriter = csv.writer(outputfile, dialect='excel')
    resultwriter.writerow(header)
    for column in result:
        resultwriter.writerow(list(column))
    outputfile.close()

inputdata = ImportCloneDataset('observeddata3.csv');
# whole course is 500 time points
# 30 time points make 2 days

#Variables:
#treatments
# treatment_cell state_killed per time pt
E_ep_killed = 0.1
E_mes_killed = 0.01
M_ep_killed = 0.01
M_mes_killed = 0.1
 

clonenum = 500
reps = 500

# fraction epithelials
# frac_ep is a vector; each row is a clone, fraction ep = E/(E+M)
EM = [2**y[1] for y in inputdata]
frac_ep_i = [y/(1+y) for y in EM]

# growth rates
# ks is a vector; each row is a clone, k is growth rate
#(modified so it matches the new time scale of 15 points per day)
ks = [y[0] for y in inputdata]     #here k is growth by days  
cell_doubling_time_i = [(15*log(2))/k for k in ks] # turn to doubling time and 15 time points / day scale

# cell numbers @ start
cells_at_start_i = [y[2] for y in inputdata]

#plasticity
#plast = 0.8
plast = 0.2

#figure out E->M and M->E for each clone. Store in vector of [(M->E,E->M)]
pr_trans_i = [((eq*plast)/(1+eq),(plast/(1+eq))) for eq in [(fep/(1-fep)) for fep in frac_ep_i]]

#set up clones; vector, one row per clone, each clone is (E cells, M cells)
clones_i = [(fep*cells_at_start_i[i],(1-fep)*cells_at_start_i[i]) for i,fep in enumerate(frac_ep_i)]


#_________________________________________________________________________________
def Treatment(clones,E_killed,M_killed,pr_trans,celldiv):
    #growth
    clones = [(cl[0]*2**(float(1)/div),cl[1]*2**(float(1)/div)) for (cl,div) in zip(clones,celldiv)]
    #death
    clones = [(cl[0]*(1-E_killed),cl[1]*(1-M_killed)) for cl in clones]
    #differentiation
    clones = [(cl[0] + cl[1]*pl[0]/div - cl[0]*pl[1]/div, cl[1] + cl[0]*pl[1]/div - cl[1]*pl[0]/div) for (pl,cl,div) in zip(pr_trans,clones,celldiv)]
    return clones

def Rest(clones,pr_trans,celldiv):
    #growth
    clones = [(cl[0]*2**(float(1)/div),cl[1]*2**(float(1)/div)) for (cl,div) in zip(clones,celldiv)]
    #differentiation
    clones = [(cl[0] + cl[1]*pl[0]/div - cl[0]*pl[1]/div, cl[1] + cl[0]*pl[1]/div - cl[1]*pl[0]/div) for (pl,cl,div) in zip(pr_trans,clones,celldiv)]
    return clones

#___________________________________________________________________________________


frac_ep_r = []
cells_at_end_r = []
fc_cells_r = []
nclones = len(clones_i)
for  re in range(reps):
    #select clones
    i_s = [random.randrange(nclones) for x in range(clonenum)] #randomly select clone IDs
    clones = [clones_i[i] for i in i_s]
    pr_trans = [pr_trans_i[i] for i in i_s]
    frac_ep = [frac_ep_i[i] for i in i_s]
    cell_doubling_time = [cell_doubling_time_i[i] for i in i_s]
    cells_treat_start = [sum(x) for x in clones]

    # run simulation
    for repeats in range(9):
        for x in range(1):
            for t in range(30):
                clones = Treatment(clones,E_ep_killed,E_mes_killed,pr_trans,cell_doubling_time)
                
            for t in range(20):
                clones = Rest(clones,pr_trans,cell_doubling_time)

        for x in range(1):
          
            for t in range(30):
                clones = Treatment(clones,M_ep_killed,M_mes_killed,pr_trans,cell_doubling_time)
                
            for t in range(20):
                clones = Rest(clones,pr_trans,cell_doubling_time)

    

    # compute the relative cell numbers (end / start)
    cells_at_end = [sum(x) for x in clones]
    fc_cells_r.append([cells[1]/cells[0] for cells in zip(cells_treat_start,cells_at_end)])
    cells_at_end_r.append(cells_at_end)

    # append frac ep to results
    frac_ep_r.append(frac_ep)

    

    

# for part (a) collect results of the simulation

#for each simulated tumor
# Bin by fraction epithelial, find the quantiles .1, .5, .9.
nbin = 10
Quant1 = []
Quant5 = []
Quant9 = []
edges = [float(x)/100 for x in range(0,100,100/nbin)]
edges.append(float(1))
bw = (edges[1]-edges[0])/2
cntrs = [round(x + bw,3) for x in edges[0:-1]]
for n in range(len(cntrs)):
    medians_temp = []
    for r in range(reps):
        fctmp = [x for i,x in enumerate(fc_cells_r[r]) if frac_ep_r[r][i] >= edges[n] and frac_ep_r[r][i] <= edges[n+1]]
        fctmp = sorted(fctmp)
        L = len(fctmp)
        if L >1: # append the median
            medians_temp.append(fctmp[int(round(.5*L))])
        elif L == 1: # or if there is just 1
            medians_temp.append(fctmp[0])
        
    medians_temp.sort()
    L = len(medians_temp)
    Quant1.append(medians_temp[int(round(.1*L))])
    Quant5.append(medians_temp[int(round(.5*L))])
    Quant9.append(medians_temp[int(round(.9*L))])

    
summary = [x for x in zip(cntrs,Quant5,Quant9,Quant1)]    



# for part (b) collect results all the simulations... figure out a representitive one later

cells_at_end = sorted([sum(c) for c in cells_at_end_r])
#compute the median, .1 and .9 quantiles of cells at end
Q1 = cells_at_end[int(round(.1*reps))]
Q5 = cells_at_end[int(round(.5*reps))]
Q9 = cells_at_end[int(round(.9*reps))]
summary.append([float('nan'), Q5, Q9, Q1])

WriteToOutput(['fracep','Median','Quart90','Quart10'],summary,'simresult16_12_19.csv')


