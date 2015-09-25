# Copyright 2015, Chen Sun (bbsunchen at outlook.com)
from sys import argv
from scipy.stats import beta

frequence_list = []

with open(argv[1]) as input_file:
    for line in input_file:
        line = line.strip()
        if line.startswith('rs#'):
            continue
        columns = line.split(' ')
        #print line
        ref_freq = float(columns[11])
        oth_freq = float(columns[14])
        if ref_freq >= 0.05 and ref_freq <= 0.95:
            frequence_list.append(ref_freq)
        if oth_freq >= 0.05 and oth_freq <= 0.95:
            frequence_list.append(oth_freq)
            
# direct fit Beta distribution
#print beta.fit(frequence_list)

# fit Beta distribution with frequency [0.05, 0.95], more precise
a,b,l,s=beta.fit(frequence_list, floc=0.04999999999999999, fscale=0.9000000000000000)
#print str(a)+'\t'+str(b)
print '{}\t{}'.format(a,b)
# fit Beta distribution with some infer
#print beta.fit(frequence_list, floc=0, fscale=1)
