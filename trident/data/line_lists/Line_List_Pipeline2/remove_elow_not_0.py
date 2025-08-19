import numpy as np
import csv
import pickle
import copy

#This file just loads the MortCash2.txt line list and removes any line with elow!=0 (ie, non-resonant lines). Because the new object is filled with 0-values, we call it, "line_dictionary0" (as opposed to 0 denoting the "first" object)


data=[]
line_dictionary0, line_dictionary={}, {}
for line in open('MortCash2.txt').readlines():
  if line.startswith('#'):
    continue
  online=line.rstrip().split()
  ion=' '.join(online[0:2])
  if ion=='':
    continue
  line_dictionary[ion], line_dictionary0[ion]=[], []    	#creating blank dictionaries
  data.append([ion]+online[2:])

for ion in line_dictionary.keys():
  for line in data:
    if ion in line:
      line_dictionary[ion].append(line)  
      
for ion, lines in line_dictionary.items():
  if ion=='':
    continue
  data0=[]
  for line in lines:
    if (line[-1]=='N/A' or float(line[-1])==0.):    #'N/A' catches only three lines in O I that we are confident should be there
      data0.append(line)
  line_dictionary0[ion]=data0




ion_column, lambda_column, gamma_column, f_column, source_column, elow_column=[], [], [], [], [], []
columns=[ion_column, lambda_column, gamma_column, f_column, source_column, elow_column]
for ion, lines in line_dictionary0.items():
  if ion=='':
    continue
  for line in lines:
    for i in range(6):
      columns[i].append(line[i])
   
  

   
datafile='MortCash0.txt'
data=np.column_stack(columns)
#data=np.column_stack([np.array(ion_column), lambda_column, gamma_column, f_column, source_column, elow_column])
np.savetxt(datafile, data, fmt='%s')










