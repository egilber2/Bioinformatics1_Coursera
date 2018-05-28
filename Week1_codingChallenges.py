# -*- coding: utf-8 -*-
"""
Created on Mon May 28 05:28:25 2018

@author: ericg
"""

#%%

    
    

#%%
# pattern count 
def readData(filename):
    with open(filename, 'r') as f:
        Text = f.readline()
        Pattern = f.readline()
        return Text.strip(), Pattern.strip()
    
from DNA_funcs import PatternCount

Text, Pattern = readData('dataset_2_7.txt')
PatternCount(Text, Pattern) #35


#%%
# Frequent Words challenge
def readData(filename):
    with open(filename, 'r') as f:
        Text = f.readline()
        k = f.readline()
        return Text.strip(), int(k.strip())

from DNA_funcs import FrequentWords

Text, k = readData('dataset_2_10.txt')

FrequentWords(Text, k)  #TAATCCGCGGTA

#%%
def readData(filename):
    with open(filename, 'r') as f:
        p = f.readline()
        return p.strip()

from DNA_funcs import revComplement

p = readData('dataset_3_2.txt')

revComplement(p)

#%%
#pattern matching
def readData(filename):
    with open(filename, 'r') as f:
        pattern = f.readline()
        genome = f.readline()
        return pattern.strip(), genome.strip()
    
from DNA_funcs import patternMatch

pattern, genome = readData('dataset_3_5.txt')

patternMatch(pattern, genome)

#%%
# computing frequencies
def readData(filename):
    with open(filename, 'r') as f:
        Text = f.readline()
        k = f.readline()
        return Text.strip(), int(k.strip())
    
from DNA_funcs import PatternToNumber, ComputingFrequencies

Text, k = readData('dataset_2994_5.txt')

ComputingFrequencies(Text, k)

#%%
#pattern to number
def readData(filename):
    with open(filename, 'r') as f:
        pattern = f.readline()
        return pattern.strip()
    f.closed

from DNA_funcs import PatternToNumber, symbolToNumber

pattern = readData('dataset_3010_2.txt')

PatternToNumber(pattern)

#%%
#Number to pattern
def readData(filename):
    with open(filename, 'r') as f:
        index = f.readline()
        k = f.readline()
        return int(index.strip()), int(k.strip())
        
from DNA_funcs import NumberToPattern

index, k = readData('dataset_3010_5.txt')

NumberToPattern(index, k)

#%%
# clump finding
def readData(filename):
    with open(filename, 'r') as f:
        genome = f.readline()
        return genome.strip()

from DNA_funcs import BetterClumpFinding, ComputingFrequencies
from DNA_funcs import PatternToNumber, NumberToPattern, symbolToNumber

genome = readData('dataset_4_5.txt')
k = 11
t = 17
L = 589

BetterClumpFinding(genome, k, t, L)



