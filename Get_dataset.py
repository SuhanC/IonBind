import pandas as pd
import random

def get_binding_site(bind_tsv,protname):
    # bindinfo = pd.read_table('/Users/suhancho/data/Uniprot_metalbinding_challenge/chebi/Zn(2+).tsv',sep='\t')
    bindinfo = pd.read_table(bind_tsv,sep='\t')
    binding_sites = bindinfo[bindinfo.Accession==protname]
    binding_sites = binding_sites.drop_duplicates(subset = ['Accession','Name'])
    return(binding_sites)

def get_dataset(test_pssm,binding_sites):
    pos_label=[]
    for i in range(len(binding_sites)):
        b_site = binding_sites.Position.tolist()[i]
        windowed = test_pssm.iloc[b_site-4 : b_site+5,:].astype(float)
        pos_label.append(windowed)

    neg_label=[]
    nb_site = list(range(len(test_pssm)))
    b_site = binding_sites.Position.tolist()
    for b in b_site:
        nb_site.remove(b)
    nb_site_sampled = random.sample(nb_site,len(pos_label))
    for nb in nb_site_sampled:
        neg_label.append(test_pssm.iloc[nb-4 : nb+5,:].astype(float))

    return(pos_label,neg_label)

    