#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pickle
import numpy as np
from random import shuffle
from tkinter import *

def RandomWeights(n):
    pre_tril = np.random.randint(1,3*n,size=(n,n))
    zero_tril = np.tril(pre_tril,-1)
    return zero_tril+zero_tril.T

def getLength(path):
    temp = 0
    for i,j in zip(path[:-1],path[1:]):
        temp += WEIGHT[i,j]
    return temp

def Populate(pop):
    for i in range(POP_SIZE):
        chromo = list([ENTRY,SINK])
        for gene in range(NET_SIZE-2):
            chromo.insert(-1,np.random.randint(NET_SIZE))
        pop.append(list(chromo))

def sort_fst(self):
        self.sort(key=getLength)
    
def SelectTopBot(pop,prcnt):
    top_notch = int(np.ceil(POP_SIZE * prcnt))
    return list(pop[:top_notch-1]),list(pop[top_notch-1:])

def merge_half(chr_a,chr_b):
    mrg_pnt = np.random.random_integers(NET_SIZE-2)
    tmp1, tmp2 = list(chr_a[:mrg_pnt]), list(chr_b[:mrg_pnt])
    tmp1.extend(chr_b[mrg_pnt:])
    tmp2.extend(chr_a[mrg_pnt:])
    if getLength(tmp1) > getLength(tmp2):
        return tmp1
    return tmp2

def mutate(chromo):
    tmp = list(chromo)
    ind = np.random.random_integers(NET_SIZE-2)
    tmp[ind] = np.random.randint(NET_SIZE)
    return tmp

def experiment(pop):
    output = list()
    for chromo in pop:
        output.append(str(chromo))
    
    sort_fst(pop)
    tmp, bot = SelectTopBot(pop,0.25) # 25-50-25 : parents-kids-weirdos
    result = list(tmp)
    for i in range(len(tmp)-1):
        for j in range(i,len(tmp)):
            result.append(merge_half(tmp[i],tmp[j]))    # splicing
    while(len(result)<POP_SIZE):
        result.append(bot.pop(0))
    
    for i in range(len(result)):
        output[i] = output[i]+' -> '+str(result[i])
        result[i] = mutate(result[i])
        output[i] = output[i]+' -> '+str(result[i])+' | '+ \
          str(getLength(pop[i]))+'->'+str(getLength(result[i]))
    sort_fst(result)
    for line in output:
        print(line)  # select, splice, mutate
    print('GEN_ENDED')
    return result

if __name__ == '__main__':
    print('START')
    NET_SIZE = 10
    POP_SIZE = NET_SIZE*2
    WEIGHT = RandomWeights(NET_SIZE)
    ENTRY = 0
    SINK = 9
    PASS_LIMIT = 30
    NUM_OF_EXPS = 100
    
    test = list()
    major_stat = list()
    for time in range(10):
        test = list()
        Populate(test)
        
        sort_fst(test)
        passes = 0
        stat = list()
        for i in range(NUM_OF_EXPS):    # SCENARIO
            pre_pop = list(test)
            stat.append([i,str(test[0]),getLength(test[0])])
            new_pop = list(experiment(test))
            if getLength(new_pop[0])<getLength(pre_pop[0]):
                test = list(new_pop)
                passes = 0
            else:
                passes += 1
                if passes >= PASS_LIMIT:
                    print('REACHED NO-CHANGE LIMIT: {}/{}'.format(i,NUM_OF_EXPS))
                    break
            #print(getLength(pre_pop[0]),'\t',getLength(test[0]))
        print(str(test[0]),str(getLength(test[0])))
        print('EXP_ENDED')
    
        '''for line in stat:
            print(line)'''
        major_stat.append(stat[-1])
    print('MAJOR STAT')
    for line in major_stat:
        print(line)
