#This is the latest incarnation of a file that attempts to wrangle the Morton line data into a useable format. We want only the lines with 
import numpy as np
import csv
import pickle
import copy

morton_dict=pickle.load(open('morton_revised_raw', 'rb'))
column_indeces=[[0, 7], [19, 28], [30, 38], [59, 68], [69, 78], [79, 88]]  #locations of each data type in the file (ion, wavelength, elow, A, gamma, f)
ion_indi, l_indi, elow_indi, A_indi, gamma_indi, f_indi=[0, 7], [19, 28], [30, 38], [59, 68], [69, 78], [79, 88]

#morton_dict should give us all the data in the file organized by multiplet for any ion with lines. Basically, we find the line number in Morton of each ion that has ground-term lines. Then, we assign all lines between ion_i and ion_i+1 to ion_i. Similarly, we find each multiplet within each ion and assign all lines in the morton file between multiplet_i and multiplet_i+1 to multiplet_i. This means our multiplets contain non-data; the last multiplet in each ion contains the heading for the next ion, etc. So, the first thing we do is throw out 'blank line's and lines that don't pass a numerical test:
    
      
mults_failed_singular, mults_failed_several=[], []
#Now we remove some non-numerical data by checking each line in each multiplet in each ion for a wavelength in the right place
numerical_dict={ion:{} for ion in morton_dict.keys()}  
for ion, mults in morton_dict.items():
  for mult, data in mults.items():
    numerical_data=[]					   #4 spaces             3               2           1          0
    blankmask=np.array([False if (line=='blank line' or line=='    ' or line=='   ' or line=='  ' or line==' ' or line=='' or 'Error' in line or '(A)' in line or 'GROUND' in line or 'GRND' in line or '=' in line) else True for line in data])   
    #here we make a last-ditch effort to remove any superfluous lines from sections with only 1 entry
    numnonblanks=np.sum(blankmask)
    if numnonblanks==0:   #in this case, the line beginning the 'multiplet' is really a misidentified line of some other sort. Then we don't want any of its data, so we move to the next mult. But
      continue           #this shouldn't really happen at this stage
    elif numnonblanks==1:  #in this case, there's only one line of data in the multiplet, and we should just test that
      line=data[blankmask][0]
      try:
        np.float64(line[l_indi[0]:l_indi[1]])
        #print('executing numerical try block on mult with 1 real entry')
        numerical_data.append(line)
      except:
        #print('numerical try block failed on mult with 1 real entry')
        mults_failed_singular.append((ion, mult, line))
    else:
      data=data[blankmask]
      meanmask=np.array([False if 'Mean' in line else True for line in data]) #we have to make this mask here, because sometimes the mean is the only info we have 
      for line in data[meanmask]:											       #for that multiplet. So we only filter it out for multiplets with 1< entries. 
        try:
          np.float64(line[l_indi[0]:l_indi[1]])  
         # print('executing numerical try block on mult with several real entries')
          numerical_data.append(line)
        except:
         # print('numerical try block failed on mult with several real entries')
          mults_failed_several.append((ion, mult, line))
    numerical_dict[ion][mult.strip()]=numerical_data    #.strip() again to make multiplet names easier
        
        
#A this point, numerical_dict should have ONLY actual spectral data, organized properly by ion/multiplet. But the multiplets include both mean and non-mean data, as well as weak lines we should do away with

#pickle.dump(numerical_dict, open('morton_revised_numerical', 'wb'))      

      
#Now we create a more organized dictionary where each entry in each multiplet is a list of numerical properties
organized_dict={ion:{} for ion in numerical_dict.keys()}
no_A_dict={ion:{} for ion in numerical_dict.keys()}
no_gamma_dict={ion:{} for ion in numerical_dict.keys()}
no_f_dict={ion:{} for ion in numerical_dict.keys()}
mask_dicts=[no_A_dict, no_gamma_dict, no_f_dict]
for ion, mults in numerical_dict.items():
  for multname, data in mults.items():
    #organized_data=[]
    #for line in data:
    #  organized_data.append([line[index[0]:index[1]].strip() for index in column_indeces[1:]]+['Morton'])   #lambda, elow, A, gamma, f, source 
    organized_data=[[line[index[0]:index[1]].strip() for index in column_indeces[1:]]+['Morton'] for line in data]
    organized_dict[ion][multname]=np.array(organized_data)
    #organized_dict[ion][multname]=np.array(organized_data)
    #Now we also make dicts of masks for where there are holes in the data
    for i in range(3):
      mask_dicts[i][ion][multname]=np.array([True if line[i+2]=='' else False for line in organized_data])


#Finally, we remove any entry that is missing both gamma and A values, as these are not in Cashman; if both are missing, the line cannot be used
filled_m_dict={ion:{} for ion in organized_dict.keys()}
for ion in filled_m_dict.keys():
  filled_m_dict[ion]={mult:[] for mult in organized_dict[ion].keys()}
for ion, mults in organized_dict.items():
  for mult, data in mults.items():
    missing_needed_data_mask=no_A_dict[ion][mult]*no_gamma_dict[ion][mult]   #True and unuseable if missing both gamma and A
    #if np.sum(missing_needed_data_mask)==0
    filled_m_dict[ion][mult]=organized_dict[ion][mult][~missing_needed_data_mask]
    
#Finally, we remove any multiplets or ions for which there is no data
filled_m_dict_iterator=copy.deepcopy(filled_m_dict)
for ion, mults in filled_m_dict_iterator.items():
  for mult, data in mults.items():
    if len(data)==0:
      del filled_m_dict[ion][mult]
      
filled_m_dict_iterator=copy.deepcopy(filled_m_dict)
for ion, mults in filled_m_dict_iterator.items():
  if len(mults)==0:
    del filled_m_dict[ion]
#filled_m_dict now has all lines with recorded gamma or A values organized by ion/multiplet. 
all_dicts, dict_names=mask_dicts+[organized_dict, filled_m_dict], ['no_A', 'no_gamma', 'no_f', 'main', 'filled']
pickle.dump({dict_names[i]:all_dicts[i] for i in range(5)}, open('morton_revised', 'wb'))  
        
        
        
