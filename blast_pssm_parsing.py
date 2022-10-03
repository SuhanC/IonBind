import pandas as pd
import os

inpath='/Users/suhancho/data/Uniprot_metalbinding_challenge/PSSM/'
files = [inpath + f for f in os.listdir(inpath)]

test = open(files[1],'r').readlines()

aa_cols=['Residue','A', ' R', ' N', ' D', ' C', ' Q', ' E', ' G', ' H', ' I', ' L', ' K', ' M', ' F', ' P', ' S', ' T', ' W', ' Y', ' V']

pssm_lst=[]
for i in range(3, len(test)):
    pssm_line = test[i].split('  ')
    pssm_line = [l.strip() for l in pssm_line if l!='']
    matrix_1 = pssm_line[0:21]
    matrix_2 = pssm_line[21:41]
    pssm_lst.append(matrix_1)

pssm = pd.DataFrame(pssm_lst)
pssm.columns = aa_cols
pssm = pssm.dropna()

pssm['Residue'] =[r.split(' ')[1] for r in pssm['Residue'].tolist()]
pssm.set_index('Residue',inplace = True)