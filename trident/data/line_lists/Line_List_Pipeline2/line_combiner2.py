import numpy as np
import csv
import pickle
import copy




#This file compares and compiles info from 3 spectral files. To do so, we first load them up and put them in a common layout
#First we load pre-processed Morton:
morton_dicts=pickle.load(open('organized_morton_data', 'rb'))
morton_dict=morton_dicts['data']
#We want to copy our values for A over to gamma values wherever we are missing them:
#for ion, data in morton_dicts['data'].items():
#  for line in data:
#    if line[2]=='':
#      line[2]==line[1]

#We already made a mask for lines without gamma or A values; however, we add some trident data here that can fill in the blanks. So, we end up making NEW no_gamma/no_A masks here
#This doesn't fix lines where we have neither gamma nor A values though; done below with turning everything into array    
#(ion, wavelength, A, gamma, f, elow)

    

#Now we load up Trident
tri_data, tri_dict=[], {}
for line in open('trident_lines.txt').readlines():
  online=line.rstrip().split()
  if line.startswith("#") or len(online) < 5:
    continue
  ion=' '.join(online[0:2])
  tri_dict[ion]=[]
  tri_data.append((ion, line[15:26], '', line[30:42], line[46:58], '', 'Trident'))
  
for ion in tri_dict.keys():
  for linedata in tri_data:
    if ion in linedata:
      tri_dict[ion].append(linedata[1:])   #lambda,  A, gamma, f, elow, source

#And now Cashman
cashman_data, cashman_dict=[], {}
for line in open('cashman_lines.txt').readlines():
  if line.startswith('#'):
    continue
  online=line.rstrip().split()
  ion=' '.join(online[1:3])
  
  cashman_dict[ion]=[]
  cashman_data.append((ion, line[90:98],  '', '', line[105:113], '', 'Cashman')) #ion, lambda, A, gamma, f, elow, source
  
for ion in cashman_dict.keys():
  for linedata in cashman_data:
    if ion in linedata:
      cashman_dict[ion].append(linedata[1:]) #lambda,  A, gamma, f, elow, source
      
#Now we turn everything into arrays. Note that we have kept everything strings so that this can be done; to do numerical comparisons, will need to turn to floats inline
#cashman_dict={ion:np.array(lines) for ion, lines in cashman_dict.items()}
#tri_dict={ion:np.array(lines) for ion, lines in tri_dict.items()}
#morton_dicts['data']={ion:np.array(lines) for ion, lines in morton_dicts['data'].items()}
#morton_dicts['no_gamma']={ion:np.array(lines) for ion, lines in morton_dicts['no_gamma'].items()}
#morton_dicts['no_A']={ion:np.array(lines) for ion, lines in morton_dicts['no_A'].items()}
#no_gamma_no_A={ion:np.array(mask)*np.array(morton_dicts['no_A'][ion]) for ion, mask in morton_dicts['no_gamma'].items()}


#FINALLY, we start comparing!
#The first thing we do is combine Morton and Cashman data. If there are matching lines in Cashman, we update the Morton data with Cashman values. If there are extra lines in Cashman, we append them. Then, if there are matching lines in Trident we can fill in missing gamma/A values for Morton and Cashman. For lines in Trident with no match in MortCash, we will simply append them.



shared_ions_mortcash=[ion for ion in cashman_dict.keys() if ion in morton_dict.keys()]
exclusive_ions_m=[ion for ion in morton_dict.keys() if ion not in cashman_dict.keys()]  #These ions are only in Morton 
exclusive_ions_cash=[ion for ion in cashman_dict.keys() if ion not in morton_dict.keys()] #These ions are only in Cashman, so we can't use due to lack of gamma/A values. But we do have one of these ions in Trident (Ar I), so maybe we'll use Trident data for that one
#shared_lines, exclusive_lines={}


d_lambda_threshold=0.2
combined_dict={ion:[] for ion in shared_ions_mortcash}
picked_lines=[]
doublepicked=[]
incashnotmort={ion:[] for ion in cashman_dict.keys()}
for ion in shared_ions_mortcash:
  #We want to replace line data present in Morton with that present in Cashman, where available. So check each Morton line-->if in Cashman, replace values.
  #If not in Cashman, keep Morton values. For values that ARE in Cashman, but not Morton, we don't have full data. We check Trident file for these values later
  #If len(data)==1, make it a list--> an array, to get array len 1 which is iterable, when array len 0 is not. We can't just make them all array([]) first, as it turns non-len1 lists into len1 arrays, for instance array([[1, 2, 3]])
  c_data, m_data=np.array(cashman_dict[ion], dtype='<U12'), np.array(morton_dict[ion], dtype='<U12')
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
    if min(deltas)<d_lambda_threshold:          #Check if close enough
      #print('1', deltas.shape)
      m_index=np.squeeze(np.where(deltas==min(deltas)))  #This tells us which Morton line it is
      try:
        len(m_index)
        print(ion, m_index)
      except:
        pass
      #replace the data for this line with lambda, f values from Cashman data
      data[m_index, 0], data[m_index, 3], data[m_index, -1]=c_data[i, 0], c_data[i, 3], 'MortCash'   
      if (ion, m_index) in picked_lines:
        doublepicked.append((ion, m_index)) #We assume all Cashman lines are either in Morton, or so far apart that they cannot be mistaken for ones that are in Morton. Because we iterate over them, each Cashman line can be matched only once. However, it is possible for the same Morton line to match two Cashman lines (doesn't occur)    
      picked_lines.append((ion, m_index))   
    else:
      #If there is no match, the line is ostensibly not in Morton-->add data from cashman to special dictionary
      incashnotmort[ion].append(c_data[i])      
  combined_dict[ion]=list(data)   #We breifly turned data into an array to do arithmetic; now turn back into list for appending   

#We now have a dictionary of shared ions between Morton/Cashman. For each shared ion, if a line was in Morton and Cashman, now has updated Cashman data. If line in Cashman but not Morton, or if we are missing gamma AND A values for Morton, do not have full data. 

#Now, add in any Morton ions not covered by Cashman
for ion in exclusive_ions_m:
  combined_dict[ion]=morton_dict[ion]
  
#for ion, lines in combined_dict.items():
#  combined_dict[ion]=np.array(lines)  
  
#Without using Trident data to fill in blanks, this is as far as we can go. Now we weed out any lines missing gamma/A values:
reduced_dict, no_gamma_dict, no_A_dict={}, {}, {}
for ion, data in combined_dict.items():
  no_A_mask, no_gamma_mask=[], []
  for line in data:
    if line[1]=='':
      no_A_mask.append(True)
    else:
      no_A_mask.append(False)
    if line[2]=='':
      no_gamma_mask.append(True)
    else:
      no_gamma_mask.append(False)
  reduced_dict[ion]=np.array(data)[~(np.array(no_gamma_mask)*np.array(no_A_mask))]   #Have to turn everything into arrays here to use mask properly
  no_gamma_dict[ion]=no_gamma_mask
  no_A_dict[ion]=no_A_mask
  
#And now we replace gamma values with A values where necessary, and add 'N/A' for any lines without elow listed
for ion, data in reduced_dict.items():
  for line in data:
    if line[2]=='':
      line[2]=line[1]
    if line[4]=='':
      line[4]='N/A'
      
pickle.dump(reduced_dict, open('mortcash_dict_new', 'wb'))

#Now we make a line list in the format Trident expects
ion_column, lambda_column, gamma_column, f_column, elow_column, source_column=[], [], [], [], [], []
for ion, data in reduced_dict.items():
  for line in data:
    ion_column.append(str(ion))
    lambda_column.append(str(line[0]))
    gamma_column.append(str(line[2]))
    f_column.append(str(line[3]))
    elow_column.append(str(line[4]))
    source_column.append(str(line[5]))

columns=[ion_column, lambda_column, gamma_column, f_column, elow_column, source_column]
datafile='MortCash.txt'
data=np.column_stack([np.array(column) for column in columns])
np.savetxt(datafile, data, fmt='%s')

#Now we add Trident data; either lines only present in Trident, or lines present in Morton/Cashman that are missing both gamma and A values. Lines that are missing gamma but not A have already had A values substituted for gamma, and the gamma values listed in Trident are not used
shared_ions_tri=[ion for ion in tri_dict.keys() if ion in combined_dict.keys()]
exclusive_ions_tri=[ion for ion in tri_dict.keys() if ion not in combined_dict.keys()]  #These ions are only in Trident
trident_picked=[]
trident_doublepicked=[]    
combined_dict_before_trident=copy.deepcopy(combined_dict)
for ion in shared_ions_tri:
  tri_lines=tri_dict[ion]
  combo_lines=combined_dict[ion]
  combo_lambda=np.float64([line[0] for line in combo_lines])
  for tri_line in tri_lines:
    tri_lambda=np.float64(tri_line[0])
    deltas2=abs(tri_lambda-combo_lambda)
    #if we find a match AND we need the A and gamma values, take gamma from Trident
    if min(deltas2)<d_lambda_threshold:
      #print('2', deltas2.shape)
      combo_index=np.where(deltas2==min(deltas2))
      #Don't use squeeze here because we want an array of len 1, not an unsized array
      #Sometimes there are two transitions that happen at the same wavelength, but have differences in other properties; these are listed seperately in
      #Morton, but not in Trident--instead, they are combined and their properties added. We prefer to keep them seperate, so if we find such an occurance
      #We simply skip it and move to the next Trident line
      if len(combo_index)>1: 
        print(ion, tri_lambda, combo_index, 'morton index longer than 1')  
        continue
      else:
        combo_index=combo_index[0]
      if (combo_lines[combo_index][1]=='' and combo_lines[combo_index][2]==''):
        #fill the blank data for this line with gamma value from Trident data
        combo_lines[combo_index][2], combo_lines[combo_index][-1]=tri_line[2], 'MortCashTri'
        if (ion, combo_index) in trident_picked:
          trident_doublepicked.append((ion, combo_index))
        trident_picked.append((ion, combo_index))
    #if there's no match, add this line to our data as a seperate line
    else:
      combo_lines.append(tri_line)      
    
  
#And now we add in Trident-exclusive ions
for ion in exclusive_ions_tri:
  combined_dict[ion]=tri_dict[ion]

#Finally, we do a last reduction, removing any lines that still lack gamma or A values (or f values! Some lines are missing f values, but I checked and this only happens to lines that are ALSO missing gamma and A values)
reduced_dict2, no_gamma_dict2, no_A_dict2={}, {}, {}
for ion, data in combined_dict.items():
  no_A_mask, no_gamma_mask=[], []
  for line in data:
    if line[1]=='':
      no_A_mask.append(True)
    else:
      no_A_mask.append(False)
    if line[2]=='':
      no_gamma_mask.append(True)
    else:
      no_gamma_mask.append(False)
  reduced_dict2[ion]=np.array(data)[~(np.array(no_gamma_mask)*np.array(no_A_mask))]   #Have to turn everything into arrays here to use mask properly
  no_gamma_dict2[ion]=no_gamma_mask
  no_A_dict2[ion]=no_A_mask

#And now we replace gamma values with A values where necessary 
for ion, data in reduced_dict2.items():
  for line in data:
    if line[2]=='':
      line[2]=line[1]
    if line[4]=='':
      line[4]='N/A'
pickle.dump(reduced_dict2, open('mortcashtri_dict_new', 'wb'))

#Now we make a line list in the format Trident expects
ion_column, lambda_column, gamma_column, f_column, elow_column, source_column=[], [], [], [], [], []
for ion, data in reduced_dict2.items():
  for line in data:
    ion_column.append(str(ion))
    lambda_column.append(str(line[0]))
    gamma_column.append(str(line[2]))
    f_column.append(str(line[3]))
    elow_column.append(str(line[4]))
    source_column.append(str(line[5]))

columns=[ion_column, lambda_column, gamma_column, f_column, elow_column, source_column]
datafile='MortCashTri.txt'
data=np.column_stack([np.array(column) for column in columns])
np.savetxt(datafile, data, fmt='%s')


#8/27/2021, email mentions errors in trident lines



"""
trident_picked=[]
trident_doublepicked=[]    
combined_dict_before_trident=copy.deepcopy(combined_dict)
mortcashtri_dict={ion:[] for ion in combined_dict}
ions_not_in_trident=('Si I', 'Cl', 'Ti', 'Co', 'Ni', 'Cu', 'Zn', 'Ga')
for ion, lines in combined_dict.items():
  if ion in (ions_not_in_trident):
    continue
  tri_lines=tri_dict[ion]
  trident_lambda=np.float64([line[0] for line in tri_lines])
  for mc_line in lines:
    mc_lambda=np.float64(mc_line[0])
    deltas=abs(mc_lambda-trident_lambda)
    #if we find a match AND we need the A and gamma values, take gamma from Trident
    if min(deltas)<d_lambda_threshold and mc_line[2]=='' and mc_line[3]=='':
      tri_index=np.squeeze(np.where(deltas==min(deltas)))
      if (ion, tri_index) in trident_picked:
        trident_doublepicked.append((ion, tri_index))
      trident_picked.append((ion, tri_index))
      mc_line[3], mc_line[-1]=tri_lines[tri_index][3], 'MortCashTri'
  mortcash_tri_dict[ion]=lines
"""

































