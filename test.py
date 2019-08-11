# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 10:33:10 2019

@author: HYEJEONG
"""
import pandas as pd
from collections import Counter
from mathematics import *
import numpy as np
import matplotlib.pyplot as plt


aminoacid = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
structure = ['h', 'e', '_']


symbollist = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
statelist = {'h', 'e', '_'}

def readfile(string):
        """
        This is a method to determine if argument is a string representing a numeric value. 
        """ 
        for kind in (str, str, str): 
            try: 
                kind(string) 
            except (TypeError, ValueError): 
                pass 
            else: 
                return True 
        else: 
            return False 

def lineseq(path): 
    """
    This method is for seperating string to word for getting symbols and states.
    """
    allstring = [] 
    with open(path) as f: 
        for line in (line.strip() for line in f):  ## line 말고 통째로 split 
         fields = line.split() 
         if fields: # non-blank line? 
             if readfile(fields[0]):
                 allstring += fields
    return allstring 

def linkerseq(path, dbtype = 1):
    if dbtype == 0: # from db file of coco
        allstring = lineseq(path)
        protein = []
        linkerstr = [] 
        for j in range(len(allstring)):      
            s = allstring[j]
            #ss = list(s)
            #print(list(s))
            #print(symbollist)
            if (set(list(s)) - set(aminoacid)) == set():
                #print('list in amino acid1')
                protein.append(list(allstring[j]))
                one = np.ones(len(allstring[j]), int)
                
                linkerstr.append(list(map(str, one)))
    elif dbtype == 1:
        allstring = lineseq(path)
        protein = []
        linkerstr = []   
   
        for j in range(len(allstring)):      
            s = allstring[j]
            if s.startswith("0") or s.startswith("x") or s.startswith("1"):
                linkerstr.append(list(allstring[j]))
                protein.append(list(allstring[j-1]))    
    return protein, linkerstr       

def proteinseq(path, dbtype = 1):
    """
    This is a method for getting a protein set. 
    output = [ [ [ protein1 ], [ protein2 ], [ protein3 ], [ protein4 ], ... ],
               [ [ structure1 ], [ structure2 ], [ structure3 ], [ structure4 ], ... ] ]
    """   
    if dbtype == 0:
        None #FASTA file reading method will be here... To be continued...      
    else:
        allstring = lineseq(path)
        protein = []
        secondstr = []   
        i = None
        for j in range(len(allstring)):    
            if allstring[j] == '<>' or allstring[j] == '>':
                protein_single = []
                secondstr_single = []
                i = j+1
                while i < len(allstring):
                    if allstring[i] == 'end' or allstring[i] == '<>' or allstring[i] == '<end>'  :
                        protein.append(protein_single)
                        secondstr.append(secondstr_single)
                        protein_single =  []
                        secondstr_single = []
                        break
                    else: #allstring[i] != '<' or allstring[i] != '>' or allstring[i] != '<>':
                        protein_single.extend(allstring[i])
                        secondstr_single.extend(allstring[i+1])
                        i += 2
    return protein, secondstr         

def getproteinset(path):
    proteinset = proteinseq(path) 
    return proteinset                        
    
def getlinkerset(path):
    linkerset = linkerseq(path) 
    return linkerset 
## Counting and mathematics from here
def count(count, item, number):
    if item not in count:
        count[item] = 0
    count[item] += number

def count1d(count, item): #Count +1 'count' if 'item' is in 'count'. 
    if item not in count:
        count[item] = 0
    count[item] += 1

def count2d(count, item1, item2): #Count +1 'count' if 'item' is in 'count', for 2-dimensional array(dictionary)
    if item1 not in count:
        count[item1] = {}
    count1d(count[item1], item2)

def norm1d(prob, item_set):  #Normalize 1-dimensional array for getting probability
    result = {} 
    prob_sum = 0.0
    
    for item in item_set:
        prob_sum += prob.get(item, 0)
    
    if prob_sum == 0:
        result[item] = 0
    else:    
        for item in item_set:
            result[item] = prob.get(item, 0) / prob_sum            
    return result
                       
def norm2d(prob, item_set1, item_set2): #Normalize 2-dimensional array for getting probability
    result = {}
    
    if prob is None:
        for item in item_set1:
            result[item] = norm1d(None, item_set2)
    for item in item_set1:
        result[item] = norm1d( prob.get(item, 0), item_set2)
    
    return result
## Finding dictionary for counting 
def initialset(proteinset):
    """
    This is a method for finding parameters by raw data sets. 
    Raw data sets are seperated as below. 
     
    """
    proteins = proteinset[0] 
    structures = proteinset[1]
    
    symbolcount = {}
    statecount = {}

    startcount = {}
    transcount = {}
    statesymbolcount = {} 
    
    for protein, structure in zip(proteins, structures):
        pre_state = {}
        count1d(startcount, structure[0])    
        for symbol, state in zip(protein, structure): #single protein
            count1d(statecount, state) 
            count1d(symbolcount, symbol)
            count2d(statesymbolcount, state, symbol) # for transmission probability           
            if pre_state != {}: 
                count2d(transcount, pre_state, state) # for emittance probability            
            pre_state = state

    statelist = list(statecount.keys())
    symbollist = list(symbolcount.keys()) 

    prob_trans = norm2d(transcount, statelist, statelist)
    prob_emit = norm2d(statesymbolcount, statelist, symbollist)
    prob_start = norm1d(startcount, statelist)
    #print(prob_start)
    return [startcount, transcount, statesymbolcount], [prob_start, prob_trans, prob_emit], \
           [statelist, symbollist]

def findLinker(proteinset):
    proteins = proteinset[0] 
    structures = proteinset[1]
    
    proteinLinker = []
    structureLinker = []
  
    
    for protein, structure in zip(proteins, structures):
        linkDistCount = {}
        var = []
        for i in range(len(structure)):
            #print(i)
            #print(structure[i:i+5])
            if structure[i:i+5] == ['1', '1', '1', '1', '1']:
                var.append(protein[i]) 
                #print(structure[i:i+5])
                #print(len(structure))
                #print(i)
                try: 
                    if structure[i+5] != '1':
                        var.append(protein[i+1])
                        var.append(protein[i+2])
                        var.append(protein[i+3])
                        var.append(protein[i+4])
                        break
                except:
                    var.append(protein[i+1])
                    var.append(protein[i+2])
                    var.append(protein[i+3])
                    var.append(protein[i+4])
                    break
            
        if var != []:
            proteinLinker.append(var)
    
    return proteinLinker

def countLinker(proteinLinker):
 #    symbollist = {'A':, 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}         
    var = {}  
    countvar = {}
    for protein in proteinLinker:
        for index, value in enumerate(protein):
            #print(value)
            
            #print(index+1)
            #print(var)
            count(var, value, (index+1)/len(protein))
            count1d(countvar, value) #count(var, value, index+1)
    
    avg = {}
    for symbol in symbollist:
        avg[symbol] = var[symbol] / countvar[symbol]
    return avg 


## main 
dbLinker = linkerseq('./dataset/TRAINING.fasta', 1) 
dbLinkerCo = linkerseq('./dataset/dna_seq_linkers_utf-co.txt', 0)

listDfl = initialset(dbLinker)[1][2]
listCo = initialset(dbLinkerCo)[1][2]

emitDfl = pd.DataFrame(listDfl)
emitCoDf = pd.DataFrame(listCo)

linkerDfl = emitDfl['1'].drop(['B', 'Z'], 0)
linkerCo = emitCoDf['1'] 

#linkerDfl['1']
print(linkerDfl) 
print(linkerCo)
'''
if set(emitCo['1'].keys()) != set(emit['1'].keys()) :
    subSet = set(emit['1'].keys()) - set(emitCo['1'].keys()) 
    #print(subSet)
    for aa in subSet:
  #      print(' emit Co !! = ', emitCo['1'][subSet] )
        emitCo['1'][aa] = 0.0
'''        
    



'''
#ax = s.groupby(symbollist).plot(kind='bar', stacked=True)
ax1 = linkerDfl.plot.bar(stacked=True, rot=0, title='Probability of amino acids found in linker (144 sets, DFLpred)')
ax1.set(xlabel='Amino acid', ylabel='Probability', ylim=[0, 0.18])
plt.savefig('freq_dfl.pdf')
plt.show()

ax2 = linkerCo.plot.bar(stacked=True, rot=0, title='Probability of amino acids found in linker (34 sets)')
ax2.set(xlabel='Amino acid', ylabel='Probability', ylim=[0, 0.18])
plt.savefig('freq_co.pdf')
plt.show()
#
#f.plot.bar(stacked=True, ='Amino acid', y='Frequency' )
#['Name'].count().unstack('Abuse/NFF').fillna(0)
'''

findlinkerDfl = findLinker(dbLinker)
findlinkerCo = findLinker(dbLinkerCo)
countLinkerDfl = countLinker(findlinkerDfl)
countLinkerCo = countLinker(findlinkerCo)
#print(countLinkerDfl)
'''
# 이건 왜? 
for var in [countLinkerDfl, countLinkerCo]:
    print(var.keys())
    if set(var.keys()) != set(symbollist) :
        print(set(var.keys()))
        subSet = set(symbollist) - set(var.keys())
        #print(subSet)
        for aa in subSet:
      #      print(' emit Co !! = ', emitCo['1'][subSet] )
            var[aa] = 0.0
'''            
           
print(countLinkerDfl)
print(countLinkerCo)

emitLinkerDfl = pd.Series(countLinkerDfl).to_frame().sort_index()
emitLinkerCoDf = pd.Series(countLinkerCo).to_frame().sort_index()
            
print(emitLinkerDfl)
print(emitLinkerCoDf)
colormark = ['gray', 'gray', \
                                    'yellow', 'yellow',\
                                    'gray',\
                                    'blue', \
                                    'red', \
                                    'gray', \
                                    'red', \
                                    'gray', 'gray',  \
                                    'orange', \
                                    'purple', \
                                    'orange', \
                                    'red', \
                                    'orange' ,'orange', \
                                    'gray', 'gray',\
                                    'orange']

a = np.array(colormark).T.tolist()
print(a)
#a= list(map(list, zip(*colormark)))
#a = list(islice(cycle(['b', 'r', 'g', 'y', 'k']), None, len(emitLinkerDfl)))
ax3 = emitLinkerDfl.plot.bar(stacked=True, rot=0, \
                             title='Location of amino acids found in linker (144 sets, DFLpred)',\
                             color='gray')

ax3.set(xlabel='Amino acid', ylim=[0.2, 0.8])
plt.savefig('freq_linker_dfl_rev1.pdf')
plt.show()

ax4 = emitLinkerCoDf.plot.bar(stacked=True, rot=0, \
                              title='Location of amino acids found in linker (34 sets)',\
                              color='gray')
ax4.set(xlabel='Amino acid',  ylim=[0.2, 0.8])
plt.savefig('freq_linker_co_rev1.pdf')
plt.show()

