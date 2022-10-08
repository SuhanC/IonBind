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




def get_dataset_padded(test_pssm,binding_sites):
    pos_label=[]
    for i in range(len(binding_sites)):
        b_site = binding_sites.Position.tolist()[i]
        pad_columns = test_pssm.columns.tolist()

        if ((b_site-4)<0 )& ((b_site+5) > len(test_pssm)) :
            front_pad_seq = ['X']*(4-b_site)
            front_pad_values = [[-0.5]*len(pad_columns)]*(4-b_site)
            front_pad = pd.DataFrame(front_pad_values)
            front_pad.index = front_pad_seq
            front_pad.columns = pad_columns
            pssm_og = test_pssm.iloc[b_site-4 : b_site+5,:].astype(float) 
            rear_pad_seq = ['X']*((5+b_site)-len(test_pssm))
            rear_pad_values = [[-0.5]*len(pad_columns)]*((5+b_site)-len(test_pssm))
            rear_pad = pd.DataFrame(rear_pad_values)
            rear_pad.index = rear_pad_seq
            rear_pad.columns = pad_columns
            windowed = pd.concat([front_pad,pssm_og,rear_pad])

        elif ((b_site-4)<0) & ((b_site+5) < len(test_pssm)):
            front_pad_seq = ['X']*(4-b_site)
            front_pad_values = [[-0.5]*len(pad_columns)]*(4-b_site)
            front_pad = pd.DataFrame(front_pad_values)
            front_pad.index = front_pad_seq
            front_pad.columns = pad_columns
            pssm_og = test_pssm.iloc[b_site-4 : b_site+5,:].astype(float) 
            windowed = pd.concat([front_pad,pssm_og])

        elif ((b_site-4)>0) & ((b_site+5) > len(test_pssm)):
            pssm_og = test_pssm.iloc[b_site-4 : b_site+5,:].astype(float) 
            rear_pad_seq = ['X']*((5+b_site)-len(test_pssm))
            rear_pad_values = [[-0.5]*len(pad_columns)]*((5+b_site)-len(test_pssm))
            rear_pad = pd.DataFrame(rear_pad_values)
            rear_pad.index = rear_pad_seq
            rear_pad.columns = pad_columns
            windowed = pd.concat([pssm_og,rear_pad])

        else:
            windowed = test_pssm.iloc[b_site-4 : b_site+5,:].astype(float)
        pos_label.append(windowed)

    ############################## Internal False labeling ###################################
    neg_label=[]
    nb_site = list(range(len(test_pssm)))
    b_site = binding_sites.Position.tolist()
    pad_columns = test_pssm.columns.tolist()

    for b in b_site:
        nb_site.remove(b)

    for b_site in nb_site : 
        if ((b_site-4)<0 )& ((b_site+5) > len(test_pssm)) :
            front_pad_seq = ['X']*(4-b_site)
            front_pad_values = [[-0.5]*len(pad_columns)]*(4-b_site)
            front_pad = pd.DataFrame(front_pad_values)
            front_pad.index = front_pad_seq
            front_pad.columns = pad_columns
            pssm_og = test_pssm.iloc[b_site-4 : b_site+5,:].astype(float) 
            rear_pad_seq = ['X']*((5+b_site)-len(test_pssm))
            rear_pad_values = [[-0.5]*len(pad_columns)]*((5+b_site)-len(test_pssm))
            rear_pad = pd.DataFrame(rear_pad_values)
            rear_pad.index = rear_pad_seq
            rear_pad.columns = pad_columns
            windowed = pd.concat([front_pad,pssm_og,rear_pad])

        elif ((b_site-4)<0) & ((b_site+5 )< len(test_pssm)):
            front_pad_seq = ['X']*(4-b_site)
            front_pad_values = [[-0.5]*len(pad_columns)]*(4-b_site)
            front_pad = pd.DataFrame(front_pad_values)
            front_pad.index = front_pad_seq
            front_pad.columns = pad_columns
            pssm_og = test_pssm.iloc[0:(b_site+5),:].astype(float) 
            windowed = pd.concat([front_pad,pssm_og])

        elif ((b_site-4)>0) & ((b_site+5) > len(test_pssm)):
            pssm_og = test_pssm.iloc[b_site-4 : b_site+5,:].astype(float) 
            rear_pad_seq = ['X']*((5+b_site)-len(test_pssm))
            rear_pad_values = [[-0.5]*len(pad_columns)]*((5+b_site)-len(test_pssm))
            rear_pad = pd.DataFrame(rear_pad_values)
            rear_pad.index = rear_pad_seq
            rear_pad.columns = pad_columns
            windowed = pd.concat([pssm_og,rear_pad])
        else:
            windowed = test_pssm.iloc[b_site-4 : b_site+5,:].astype(float)
        
        neg_label.append(windowed)

    return(pos_label,neg_label)

    