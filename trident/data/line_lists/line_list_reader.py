#The spectral line data is hard to read. This short script puts it into a dictionary, so long as any blank lines in the spectral file are blocked out with '#'
data=[]
line_dictionary={}
for line in open('MortCashTri.txt').readlines():
  if line.startswith('#'):
    continue
  online=line.rstrip().split()
  ion=' '.join(online[0:2])
  line_dictionary[ion]=[]    	#creating blank dictionary
  data.append(np.array([ion]+online[2:]))

for ion in line_dictionary.keys():
  for line in data:
    if ion in line:
      line_dictionary[ion].append(line)  
  line_dictionary[ion]=np.float64(line_dictionary[ion])

