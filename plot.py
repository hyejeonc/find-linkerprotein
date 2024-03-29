# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 18:46:38 2019

@author: HYEJEONG
"""

from mathematics import *
from statistics import *
from sequence import linkerseq

#import statistics as stat

from math import log

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

symbollist = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
statelist = ['1', 'x', '0']

dbLinker = linkerseq('./dataset/TRAINING.fasta', 1) 
dbLinkerCo = linkerseq('./dataset/dna_seq_linkers_utf-co.txt', 0) 


emit = initialset(dbLinker)[1][2]
emit_co = initialset(dbLinkerCo)[1][2]


'''
simple_emit = {'c': {'G': 0.11572760437778679, 'V': 0.05249290636400487, 'T': 0.06323469801378193, 'P': 0.0671868666396433, 'M': 0.011045804620997163, 'D': 0.07042967166599108, 'Y': 0.03090798540737738, 'N': 0.05938386704499392, 'E': 0.042359140656668015, 'Q': 0.03394811511957844, 'I': 0.030603972436157277, 'K': 0.06668017835427645, 'S': 0.0922172679367653, 'F': 0.03111066072152412, 'L': 0.06232265910012161, 'W': 0.011451155249290636, 'C': 0.023307661126874747, 'A': 0.0784353465747872, 'R': 0.03232671260640454, 'H': 0.024827725982975273}, 'e': {'G': 0.05803080308030803, 'V': 0.12788778877887788, 'T': 0.08498349834983498, 'P': 0.016226622662266228, 'M': 0.018976897689768978, 'D': 0.026402640264026403, 'Y': 0.05473047304730473, 'N': 0.030803080308030802, 'E': 0.029152915291529153, 'Q': 0.03547854785478548, 'I': 0.08085808580858085, 'K': 0.04235423542354235, 'S': 0.07288228822882288, 'F': 0.050605060506050605, 'L': 0.10451045104510451, 'W': 0.021727172717271728, 'C': 0.028052805280528052, 'A': 0.0627062706270627, 'R': 0.034653465346534656, 'H': 0.018976897689768978}, 'h': {'G': 0.04433818735057596, 'V': 0.07846120408606824, 'T': 0.04716366007389698, 'P': 0.021299717452727667, 'M': 0.021951749619647902, 'D': 0.04803303629645729, 'Y': 0.02825472723321017, 'N': 0.037383177570093455, 'E': 0.07346229080634645, 'Q': 0.042599434905455334, 'I': 0.04564225168441643, 'K': 0.08672027820039122, 'S': 0.048902412519017606, 'F': 0.04368615518365573, 'L': 0.11041078026515974, 'W': 0.01630080417300587, 'C': 0.018691588785046728, 'A': 0.12497283199304499, 'R': 0.03629645729189307, 'H': 0.025429254509889154}}
xsymbol = np.array(len(symbollist))

freq_emit = {'c': {'G': 1142, 'V': 518, 'T': 624, 'P': 663, 'M': 109, 'D': 695, 'Y': 305, 'N': 586, 'E': 418, 'Q': 335, 'I': 302, 'K': 658, 'L': 615, 'S': 910, 'C': 230, 'A': 774, 'R': 319, 'F': 307, 'W': 113, 'H': 245}, 'e': {'V': 465, 'T': 309, 'S': 265, 'F': 184, 'N': 112, 'D': 96, 'W': 79, 'G': 211, 'Q': 129, 'L': 380, 'A': 228, 'R': 126, 'Y': 199, 'I': 294, 'E': 106, 'M': 69, 'H': 69, 'K': 154, 'C': 102, 'P': 59}, 'h': {'A': 575, 'F': 201, 'D': 221, 'Q': 196, 'V': 361, 'S': 225, 'L': 508, 'K': 399, 'M': 101, 'E': 338, 'G': 204, 'T': 217, 'W': 75, 'N': 172, 'I': 210, 'C': 86, 'P': 98, 'R': 167, 'H': 117, 'Y': 130}}
#freq_trans = {'c': {'c': 8584, 'e': 729, 'h': 444}, 'e': {'e': 2906, 'c': 715, 'h': 15}, 'h': {'h': 4142, 'c': 458, 'e': 1}}


f = pd.DataFrame(freq_emit)
s = pd.DataFrame(simple_emit)
print(s)
#ax = s.groupby(symbollist).plot(kind='bar', stacked=True)
ax1 = s.plot.bar(stacked=True, rot=0, title='Simple HMM')
ax1.set(xlabel='Amino acid', ylabel='Probability', ylim=[0, 0.4])

plt.show()
#
#f.plot.bar(stacked=True, ='Amino acid', y='Frequency' )
#['Name'].count().unstack('Abuse/NFF').fillna(0)

'''
'''
for state in statelist:
    h = []
    e = []
    u = []
    for symbol in symbollist: 
        if state == 'h':
            p1 = plt.bar(symbol, simple_emit[state][symbol], label=state )
            h.append(simple_emit[state][symbol])
        elif state == 'e':
            p2 = plt.bar(symbol, simple_emit[state][symbol], label=state, bottom = h  )   
            e.append(simple_emit[state][symbol])
        elif state == 'c':
            p3 = plt.bar(symbol, simple_emit[state][symbol], label=state, bottom = e )  
            u.append(simple_emit[state][symbol])
plt.legend(p1='h', p2='e', p3='c')        
plt.show()

'''