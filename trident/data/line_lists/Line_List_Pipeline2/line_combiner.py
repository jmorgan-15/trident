import numpy as np
import csv
import pickle
import copy




#This file compares and compiles info from 3 spectral files. To do so, we first load them up and put them in a common layout
#First we load pre-processed Morton:
morton_dicts=pickle.load(open('morton_revised', 'rb'))
m_dict, filled_m_dict, no_A_dict, no_gamma_dict, no_f_dict=morton_dicts['main'], morton_dicts['filled'], morton_dicts['no_A'], morton_dicts['no_gamma'], morton_dicts['no_f']
column_indeces=[[0, 7], [19, 28], [30, 38], [59, 68], [69, 78], [79, 88]]  #locations of each data type in the file (ion, wavelength, elow, A, gamma, f)
ion_indi, l_indi, elow_indi, A_indi, gamma_indi, f_indi=[0, 7], [19, 28], [30, 38], [59, 68], [69, 78], [79, 88]
#The first thing we do is remove isotopes from the morton data:
del filled_m_dict['D I']
filled_m_dict_iterator=copy.deepcopy(m_dict)
for ion in filled_m_dict_iterator.keys():
  try:
    np.int32(ion[0])
    del filled_m_dict[ion]
  except:
    pass


    
#Now we create a version of the Morton dictionary that includes ONLY the elow=0 lines and strongest non elow=0 line for each multiplet
slim_m_dict={ion:{} for ion in filled_m_dict.keys()}
for ion in slim_m_dict.keys():
  slim_m_dict[ion]={mult:[] for mult in filled_m_dict[ion].keys()} 
for ion, mults in filled_m_dict.items():
  for mult, data in mults.items():
    elow, f=np.float64(data[:, 1]), np.float64(data[:, 4])
    if len(data)==1:
      slim_m_dict[ion][mult]=data    #If there's only one entry in this multiplet, take it. morton_prep2_2 sets this up to work correctly for multiplets with only the mean for an entry
    else:
    #Otherwise, we want all lines with elow=0   
    #3 possibilities: (yes elow=0, yes elow!=0), (yes elow=0, no elow!=0), (no elow=0, yes elow!=0)
      elow0_mask=(elow==0.) 
      elow0_data, elow_f, other_data, other_f, other_elow=data[elow0_mask], f[elow0_mask], data[~elow0_mask], f[~elow0_mask], elow[~elow0_mask]
      
      if len(other_data)==0:    #if there's only elow=0 data, just take that
        slim_m_dict[ion][mult]=elow0_data
        
         
      else:                      
        max_f_index=np.where(other_f==max(other_f))[0]            #first, find strong lines. Index is to remove list of applicable indeces from array so we can check length
        if len(max_f_index)==1:                                   #if there is a clear strongest line, just add it to the elow=0 data. This also works if there's only 1 non elow=0 line; it will just 
          slim_m_dict[ion][mult]=np.concatenate((elow0_data, other_data[max_f_index]))#return index 0 for the max and add the single needed line
          #if there are several contenders for strongest line (several lines with same f), choose one with lowest elow
        else:												
          min_elow_index=np.where(other_elow==min(other_elow[max_f_index]))     #other_elow[max_f_index] gives us the elow energies for the strongest lines; min_elow_index is the index of other_data
          slim_m_dict[ion][mult]=np.concatenate((elow0_data, other_data[min_elow_index]))  #with the lowest elow energy of the strongest line
          
      
#Now we flatten our dictionary, so only elow=0 and the strongest remaining line
flat_slim_dict={ion:{} for ion in slim_m_dict.keys()}
for ion, mults in slim_m_dict.items():
  total_data=[]
  for mult, data in mults.items():
    total_data.extend(data)
  flat_slim_dict[ion]=np.array(total_data)

#Now we load up Trident
tri_data, tri_dict0=[], {}
for line in open('trident_lines.txt').readlines():
  online=line.rstrip().split()
  if line.startswith("#") or len(online) < 5:
    continue
  ion=' '.join(online[0:2])
  tri_dict0[ion]=[]
  tri_data.append((ion, line[15:26], 'None in Trident', 'None in Trident', line[30:42], line[46:58], 'Trident'))
  
for ion in tri_dict0.keys():
  for linedata in tri_data:
    if ion in linedata:
      tri_dict0[ion].append(linedata[1:])   #lambda, elow, A, gamma, f, source

#And now Cashman
cashman_data, cashman_dict0=[], {}
for line in open('cashman_lines.txt').readlines():
  if line.startswith('#'):
    continue
  online=line.rstrip().split()
  ion=' '.join(online[1:3])
  
  cashman_dict0[ion]=[]
  cashman_data.append((ion, line[90:98], 'None in Cashman', 'None in Cashman', 'None in Cashman', line[105:113], 'Cashman'))
  
for ion in cashman_dict0.keys():
  for linedata in cashman_data:
    if ion in linedata:
      cashman_dict0[ion].append(linedata[1:]) #lambda, elow, A, gamma, f, source
      
#Now we turn everything into arrays. Note that we have kept everything strings so that this can be done; to do numerical comparisons, will need to turn to floats inline
cashman_dict={ion:np.array(lines) for ion, lines in cashman_dict0.items()}
tri_dict={ion:np.array(lines) for ion, lines in tri_dict0.items()}




#FINALLY, we start comparing!
shared_ions=[ion for ion in cashman_dict.keys() if ion in flat_slim_dict.keys()]
exclusive_ions_m=[ion for ion in slim_m_dict.keys() if ion not in cashman_dict.keys()]  #These ions are only in Morton 
exclusive_ions_cash=[ion for ion in cashman_dict.keys() if ion not in flat_slim_dict.keys()] #These ions are only in Cashman, so we can't use due to lack of gamma/A values. But we do have one of these ions in Trident (Ar I), so maybe we'll use Trident data for that one
#shared_lines, exclusive_lines={}


#First we deal with shared ions between Morton/Cashman
d_lambda_threshold=0.2
combined_dict={ion:[] for ion in shared_ions}
picked_lines=[]
doublepicked=[]
incashnotmort=[]
for ion in shared_ions:
  c_data, m_data=cashman_dict[ion], flat_slim_dict[ion]
  #We want to replace line data present in Morton with that present in Cashman, where available. So check each Morton line-->if in Cashman, replace values.
  #If not in Cashman, keep Morton values. For values that ARE in Cashman, but not Morton, we don't have full data (unless data is in Trident....)
  min_delts=[]
    #If len(data)==1, make it a list--> an array, to get array len 1 which is iterable, when array len 0 is not. We can't just make them all array([]) first, as it turns non-len1 lists into len1 arrays, for instance array([[1, 2, 3]])
  if (len(c_data)==1 and len(m_data)==1):
    lambda_m, lambda_c=np.float64([m_data[:, 0]]), np.float64([c_data[:, 0]])
    #print(ion, 'both 1')
  elif (len(c_data)==1):
    lambda_m, lambda_c=np.float64(m_data[:, 0]), np.float64([c_data[:, 0]])
    #print(ion, 'c_data 1')
  elif (len(m_data)==1):
    lambda_m, lambda_c=np.float64([m_data[:, 0]]), np.float64(c_data[:, 0])
    #print(ion, 'm_data 1')
  else:
    lambda_m, lambda_c=np.float64(m_data[:, 0]), np.float64(c_data[:, 0])
    #print(ion, 'none 1')
  data=copy.deepcopy(m_data)    #This is not strictly necessary, could just use m_data below, but it makes more sense to create a neutral data object which we fill with a mix of data
  for i in range(len(c_data)):         #for each entry in this ion in Cashman, find nearest Morton line
    deltas=abs(lambda_c[i]-lambda_m)  #Should tell us the closest Morton to this Cashman
    d_lambda=min(deltas)
    if d_lambda<d_lambda_threshold:          #Check if close enough
      m_index=np.where(deltas==d_lambda)  #This tells us which Morton line it is
      data[m_index, 0], data[m_index, 4], data[m_index, -1]=c_data[i, 0], c_data[i, 4], 'MortCash'   #replace the data for this line with Cashman data
      if (ion, m_index) in picked_lines:
        doublepicked.append((ion, m_index)) #We assume all Cashman lines are either in Morton, or so far apart that they cannot be mistaken for ones that are in Morton. Because we iterate 
      picked_lines.append((ion, m_index))   #over them, each Cashman line can be matched only once. However, it is possible for the same Morton line to match two Cashman lines (doesn't occur)    
    else: 
      incashnotmort.append((ion, lambda_c[i]))   #If there is no match, the line is ostensibly not in Morton
  combined_dict[ion]=data
    
#We now have a dictionary of shared ions between Morton/Cashman. For each shared ion, if a line was in Morton and Cashman, now has updated Cashman data. If line in Cashman but not Morton, do not have full data. 

#Now, add in any Morton ions not covered by Cashman
for ion in exclusive_ions_m:
  combined_dict[ion]=flat_slim_dict[ion]
  
for ion, lines in combined_dict.items():
  combined_dict[ion]=np.array(lines)  
#Without using Trident data to fill in blanks, this is as far as we can go. We have a spectral file!
pickle.dump(combined_dict, open('mortcash_dict2', 'wb'))


#This section for adding in Trident data. We assume that if a Cashman line is next to a Trident line, we use Trident data
incashnotmort_dict={ion:[] for (ion, lambda_) in incashnotmort} 
for ion in incashnotmort_dict.keys():
  for entry in incashnotmort:
    if ion in entry:
      incashnotmort_dict[ion].append(np.float64(entry[1]))

trident_picked=[]
trident_doublepicked=[]    
combined_dict_copy=copy.deepcopy(combined_dict)
for ion, lines in incashnotmort_dict.items():
  if ion=='Si I' or ion=='S I' or 'Cl' in ion or 'Ti' in ion or 'Co' in ion or 'Ni' in ion or 'Cu' in ion or 'Zn' in ion or 'Ga' in ion:   #No Si I data in trident, skip it
    continue
  for lambda_ in lines:
    deltas=abs(lambda_-np.float64(tri_dict[ion][:, 0]))
    if min(deltas)<d_lambda_threshold:
      tri_index=np.where(deltas==min(deltas))
      if (ion, tri_index) in trident_picked:
        trident_doublepicked.append((ion, tri_index))
      trident_picked.append((ion, tri_index))
      newdata=np.array([lambda_, 'N/A', 'N/A', tri_dict[ion][tri_index][0][3], tri_dict[ion][tri_index][0][4], 'MortCashTri']).reshape(1, 6)
      combined_dict[ion]=np.concatenate((combined_dict[ion], newdata))
      
      #np.append(combined_dict[ion], np.array([lambda_, 'N/A', 'N/A', tri_dict[ion][tri_index][0][3], tri_dict[ion][tri_index][0][4], 'MortCashTri']))
#This basically just gets us 3 oxygen lines: 936, 929, 925
      

      
ion_column, lambda_column, gamma_column, f_column, source_column, elow_column=[], [], [], [], [], []
for ion, lines in combined_dict.items():
  if ion=='':
    continue
  for line in lines:
    ion_column.append(str(ion))
    lambda_column.append(str(line[0]))
    if line[3]=='':
      gamma_column.append(str(line[2]))
    else:
      gamma_column.append(str(line[3]))
    f_column.append(str(line[4]))
    source_column.append(str(line[5]))
    elow_column.append(str(line[1]))
  

   
  
datafile='MortCash2.txt'
data=np.column_stack([np.array(ion_column), lambda_column, gamma_column, f_column, source_column, elow_column])
np.savetxt(datafile, data, fmt='%s')

"""
#to easily read data from this file:
#data=[]
#line_dictionary={}
#for line in open('MortCash2.txt').readlines():
#  if line.startswith('#'):
#    continue
#  online=line.rstrip().split()
#  ion=' '.join(online[0:2])
#  line_dictionary[ion]=[]    	#creating blank dictionary
#  data.append([ion]+online[2:])
#
#for ion in line_dictionary.keys():
#  for line in data:
#    if ion in line:
#      line_dictionary[ion].append(line) 
"""













