#This is the latest incarnation of a file that attempts to wrangle the Morton line data into a useable format. We want only the lines with 
import numpy as np
import csv
import pickle

morton_data0=list(csv.reader(open('morton_lines.txt', 'r'), delimiter="\t"))
morton_data=np.array(['blank line' if ((len(entry)==0) or (len(entry[0])==1)) else entry[0] for entry in morton_data0[12:]]) #skip header of file and make blank lines explicit
column_indeces=[[0, 7], [19, 28], [59, 68], [69, 78], [79, 88]]  #locations of each data type in the file (ion, wavelength, A, gamma, f)
ion_indi, l_indi, A_indi, gamma_indi, f_indi=[0, 7], [19, 28], [59, 68], [69, 78], [79, 88]

#First, lets find the line numbers where each ion state begins. This part leaves us with different isotopes as well, which we may not want. But it's easier to produce data for each isotope as though they were seperate elements, and then remove isotope entries later
found_ions=[]
ion_indeces=[]
file_indeces=[]
for i in range(len(morton_data)):
  line=morton_data[i]
  if (line!='blank line'):
      if ((('GROUND' in line) or ('GRND' in line)) and ('No ground' not in line) and (line[:7] not in found_ions)):
      #if ((('GROUND' in line) or ('GRND' in line) or ('No ground' in line)) and (line[:7] not in found_ions)):
      #The lines with ion names in them also have 'GROUND' or 'GRND' to specify things about the ground state. 'No ground' means the ion has no UV ionic transitions from ground state
        found_ions.append(line[:7])
        ion_indeces.append(i)
        file_indeces.append(i+13)   #this tells us what line of the file the ion is on





#Now, lets turn each ion/isotope into a dictionary key, with the lines between each as values. Luckily each ion's data is grouped together in the file
raw_dict={}
for i in range(len(found_ions)-1):
  start_index, end_index=ion_indeces[i], ion_indeces[i+1]
  raw_dict[found_ions[i]]=[morton_data[start_index:end_index], [start_index+13, end_index+13]] #each entry in the dictionary is all lines between the ion/isotope, then the lines in the file where the data is

#Now we need to group the data together in multiplets. To do this, we need to find where each multiplet starts; not as easy as you'd think, because some "multiplets" have only 1 line, and so are formatted differently. We avoid removing anything from our "data" until the end. 
found_mults, mult_indeces=[], []
for i in range(len(morton_data)):
  line=morton_data[i]
  #if ((line!='blank line') and (i not in ion_indeces) and ('No ground' not in line)):
  if ((line!='blank line') and (i not in ion_indeces) and ('No ground' not in line) and ('GROUND' not in line) and ('GRND' not in line)):
  #Sometimes the ion is repeated before continuing with the multiplets; when this happens, the ground state of the ion is also repeated, and in the same column(s) as when the state for multiplets are used. This leads to those headings being misidentified as multiplet headings. By eschewing those lines that say something about the ground state, we can avoid misinterpreting repeated ion headings with new multiplet headings
  #if ((line!='blank line') and (i not in ion_indeces)):
    try:
      z=int(line[7])
      found_mults.append(line[7:27])
      mult_indeces.append(i)      #add 13 to get Morton file line number
      print(z, line[7:27])
    except:
      pass
      
found_mults, mult_indeces=np.array(found_mults), np.array(mult_indeces)


#Now mult_indeces holds the line number in the Morton file of all lines that start with an atomic state; these lines are followed either by a multiplet mean or a single transition line. For multiplets we want to take any entries with 0 lower energy, and then the highest f-value non-0. To do this, we construct a dictionary of dictionaries; ion-->multiplet-->value(s)
mult_dict0={ion:{} for ion in raw_dict.keys()}
for i in range(len(found_ions)-1):
  current_ion, next_ion=found_ions[i], found_ions[i+1]
  current_i_index, next_i_index=ion_indeces[i], ion_indeces[i+1]   #add 13 to get Morton file line number
  if current_ion=='Ni I   ':
    next_i_index=6990-13					     #Ni I is the only entry in the file with an, "excited lower term' section. Not yet sure how to handle it, so we hardcode in a skip
  ion_mult_mask=(current_i_index<mult_indeces) & (mult_indeces<next_i_index)   #This should return only the indeces of the multiplets belonging to this ion
  for j in range(len(mult_indeces[ion_mult_mask])-1):
    mult_data=[]   #Append each line in the multiplet to this list
    current_mult, next_mult=found_mults[ion_mult_mask][j], found_mults[ion_mult_mask][j+1]
    current_m_index, next_m_index=mult_indeces[ion_mult_mask][j], mult_indeces[ion_mult_mask][j+1]
    mult_dict0[current_ion][current_mult]=morton_data[current_m_index+1:next_m_index] #This leaves off the very first line (the one that has the beginning of the current multiplet) and the very last line (the one that has the beginning of the next multiplet)
  mult_dict0[current_ion][next_mult]=morton_data[next_m_index+1:next_i_index]   #This grabs the last multiplet in each ion, which ends not with the beginning of another multiplet, but the beginning of another ion

#At this point, mult_dict should be a dictionary with each ion/multiplet combo as keys, and the lines of the Morton file between between each as values. This means the multiplets will contain non-data, like sections for ions that do not have any lines of interest to us (if ion2 has no ground term lines, 'ion 2 No ground term lines' will appear as an entry in ion1s last multiplet). The last mutliplet of each ion also contains the beginning of the section for the next ion. But both of these should be easy to detect by looking for numerical data

#Just making keys easier:
mult_dict={ion.strip():mults for ion, mults in mult_dict0.items()}
pickle.dump(mult_dict, open('morton_revised_raw', 'wb'))
    
    
      
  



"""
#Now we remove some non-numerical data:
numerical_dict={}  
for ion, data_and_index in raw_dict.items():
  numerical_data=[]
  raw_data, start_index=data_and_index[0], data_and_index[1][0]
  for i in range(len(raw_data)):
    try:
      np.float64(raw_data[i][l_index[0]:l_index[1]])   #for each line in the file currently assigned to this ion, check if it has a numerical value where the wavelength column should be. If it does, append it to numerical list
    except:
      continue
    numerical_data.append([raw_data[i], start_index+i])   #this saves the numerical data and the line of the file it is on
  numerical_dict[ion]=numerical_data
#each entry in numerical_dict is an ion. The value for each is list of length-2 lists: a string representing the line from the file that passed the wavelength test, and the line number in the Morton file of that line


#Now we check each column of each line for correct numerical data. We also construct dictionaries of masks: l_mask reflects where there is or is not numerical data in the wavelength column (this should be only true at this point in the script).
l_mask_dict, A_mask_dict, gamma_mask_dict, f_mask_dict={}, {}, {}, {}
mask_dicts=[l_mask_dict, A_mask_dict, gamma_mask_dict, f_mask_dict]
for ion, num_data_index in numerical_dict.items():
  for j in range(len(column_indeces-1)):  #don't care about ion column (or wavelenth, hopefully, but leaving for tetsing purposes)
"""    
      
        
      
        
        
        
        
