#This is the latest incarnation of a file that attempts to wrangle the Morton line data into a useable format. We want only the lines with elow=0. We also create masks for where there are holes in the data (gamma, A values)
import numpy as np
import csv
import pickle
import copy

morton_data0=list(csv.reader(open('morton_lines.txt', 'r'), delimiter="\t"))
morton_data=np.array(['blank line' if ((len(entry)==0) or (len(entry[0])==1)) else entry[0] for entry in morton_data0[12:]]) #skip header of file and make blank lines explicit
column_indeces=[[0, 7], [19, 28], [59, 68], [69, 78], [79, 88], [30, 38]]  #locations of each data type in the file (ion, wavelength, A, gamma, f, elow)
ion_indi, l_indi, A_indi, gamma_indi, f_indi, elow_indi=[0, 7], [19, 28], [59, 68], [69, 78], [79, 88], [30, 38]

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

#Previously we tried to group lines into multiplets within ions, so that we could select the striongest elow!=0 line from each. But we eventually decided to keep only elow=0 lines. Also, the multiplet part of the code led to weird errors where many ions inherited the last listed multiplet from the ion above them. So we're totally nixxing that part and going for something much simpler; just remove extreneous lines, multiplet means, and clean up the data already in raw_dict.


numerical_dict={}
for ion, lines_index in raw_dict.items():
  lines=lines_index[0]                                       #4 spaces        3             2             1           0
  blankmask=np.array([False if (line=='blank line' or line=='    ' or line=='   ' or line=='  ' or line==' ' or line=='' or 'Error' in line or '(A)' in line or 'GROUND' in line or 'GRND' in line or '=' in line) else True for line in lines]) 
  #These are the various ways in which blank lines can be represented, followed by substrings found in headers etc
  #We're going to try to assemble these masks all at once. Then we can combine them and only the lines that make it through every mask are added to our numerical dictionary
  numerical_mask=[]
  
  for i in range(len(lines)):
    line=lines[i]
     
    #this block checks for wavelength value and elow in the right spot and removes most non-numerical data   
    try:
      lambda_, elow=float(line[l_indi[0]:l_indi[1]]), float(line[elow_indi[0]:elow_indi[1]])  
    except:
      numerical_mask.append(False)
      continue
            
    #this block checks if elow is 0  
    if elow!=0.:
      numerical_mask.append(False)
      continue
          
    #this block checks if this line is a multiplet mean  
    if 'Mean' in line: 
      #There are two kinds of multiplet means with elow=0: ones where each line in multiplet has elow=0, in which case we do not want the mean. The other kind lists only the mean for the whole
      #multiplet, so we DO take it. We check this by seeing if the multiplet mean line is followed by an individual line or by the next multiplet/a blank line
      if elow==0.:
        nextline=lines[i+1]   
        try:
          lambda_next=float(nextline[l_indi[0]:l_indi[1]])
          numerical_mask.append(False)
          continue
        except:
          numerical_mask.append(True)
          continue
      else:
        mult_mean_mask.append(False)
        continue
    #If the line has lambda and elow in the right spot, elow==0, and 'Mean' isn't in line, then it is a valid peice of data. To have made it past all the continues above it should be a real line
    numerical_mask.append(True)
  print(len(blankmask), len(numerical_mask))      
  final_mask=blankmask*numerical_mask   
  numerical_dict[ion.strip()]=lines[final_mask]
  
#Now we want to remove any isotopes from the Morton data, and turn our line-strings into arrays of data. We also make masks of where there are holes in the data
del numerical_dict['D I']
ions=list(numerical_dict.keys())
for ion in ions:
  try:
    np.int32(ion[0])
    del numerical_dict[ion]
  except:
    pass    
    
organized_dict, no_A_mask_dict, no_gamma_mask_dict={}, {}, {}
for ion, data in numerical_dict.items():
  no_A_mask, no_gamma_mask=[], []
  organized_data=[[line[index[0]:index[1]].strip() for index in column_indeces[1:]]+['Morton'] for line in data]
  for line in organized_data:
    if line[2]=='':
      no_A_mask.append(True)
    else:
      no_A_mask.append(False)
    if line[3]=='':
      no_gamma_mask.append(True)
    else:
      no_gamma_mask.append(False)
  no_A_mask_dict[ion]=no_A_mask
  no_gamma_mask_dict[ion]=no_gamma_mask
  organized_dict[ion]=organized_data


#lastly, we remove any ions for which we have NO elow=0 lines in Morton. These ions may be handled exclusively by Cashman/Trident later
organized_dict_iterator=copy.deepcopy(organized_dict)
for key, value in organized_dict_iterator.items():
  if value==[]:
    del organized_dict[key]
    del no_A_mask_dict[key]
    del no_gamma_mask_dict[key]
morton_dict={'data':organized_dict, 'no_A':no_A_mask_dict, 'no_gamma':no_gamma_mask_dict}
pickle.dump(morton_dict, open('organized_morton_data', 'wb'))



"""  
#This version uses several masks and then combines them, but then it has issues redefining variables. Easier to do one mask from the start like above and not append 'True' to the maks until every test is passed
numerical_dict={}
for ion, lines_index in raw_dict.items():
  lines=lines_index[0]                                       #4 spaces        3             2             1           0
  blankmask=np.array([False if (line=='blank line' or line=='    ' or line=='   ' or line=='  ' or line==' ' or line=='' or 'Error' in line or '(A)' in line or 'GROUND' in line or 'GRND' in line or '=' in line) else True for line in lines]) 
  #These are the various ways in which blank lines can be represented, followed by substrings found in headers etc
  #We're going to try to assemble these masks all at once. Then we can combine them and only the lines that make it through every mask are added to our numerical dictionary
  elow0_mask, mult_mean_mask, lambda_mask=[], [], []
  
  for i in range(len(lines)):
    line=lines[i]
     
    #this block checks for wavelength value in the right spot and removes most non-numerical data   
    try:
      lambda_, elow=float(line[l_indi[0]:l_indi[1]]), float(line[elow_indi[0]:elow_indi[1]])  
      lambda_mask.append(True)
    except:
      lambda_mask.append(False)
            
    #this block checks if elow is 0  
    if elow==0.:
      elow0_mask.append(True)
    else:
      elow0_mask.append(False)
          
    #this block checks if this line is a multiplet mean  
    if 'Mean' in line: 
      #There are two kinds of multiplet means with elow=0: ones where each line in multiplet has elow=0, in which case we do not want the mean. The other kind lists only the mean for the whole
      #multiplet, so we DO take it. We check this by seeing if the multiplet mean line is followed by an individual line or by the next multiplet/a blank line
      if elow==0.:
        nextline=lines[i+1]   
        try:
          lambda_next=float(nextline[l_indi[0]:l_indi[1]])
          mult_mean_mask.append(False)
        except:
          mult_mean_mask.append(True)
      else:
        mult_mean_mask.append(False) 
  print(len(blankmask), len(elow0_mask), len(mult_mean_mask), len(lambda_mask))      
  numerical_mask=blankmask*elow0_mask*mult_mean_mask*lambda_mask     
  numerical_dict[ion.strip()]=lines[numerical_mask]
""" 
    
  




      
        
      
        
        
        
        
