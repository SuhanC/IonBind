import pandas as pd

def getdata():
    datpath = '/Users/suhancho/data/Uniprot_metalbinding_challenge/'
    train_pos = pd.read_table(datpath+'POS_TRAIN_FULL.tsv')
    chebi = pd.read_table(datpath+'ChEBI-IDs_for_metal_binding.tsv')
    return(train_pos,chebi)
train_pos,chebi = getdata()
train_pos = pd.merge(train_pos,chebi,on = 'ChEBI-ID')

OUTPATH='/Users/suhancho/data/Uniprot_metalbinding_challenge/chebi/'
for cb in train_pos['ChEBI-ID'].unique().tolist():
    train_tmp = train_pos[train_pos['ChEBI-ID']==cb]
    train_tmp.to_csv(OUTPATH+train_tmp['Name'].values[0].replace(' ','')+'.tsv',sep='\t',index = False)