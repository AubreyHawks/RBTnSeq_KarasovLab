import pickle
import csv

###############################################
#path to the _pos_to_gene.cpk dictionary
###############################################
cpk = 'demo_pos_to_gene.cpk'

with open(cpk, 'rb') as pickles:
    cpk_dict = pickle.load(pickles)

##############################################
#rename output file here
##############################################
with open('cpk_dict.csv', 'w') as dict_csv:
    for key in cpk_dict.keys():
        dict_csv.write("%s,%s\n"%(key,cpk_dict[key]))
