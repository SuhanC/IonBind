from torch.utils.data import DataLoader, Dataset
import torch.nn.functional as F
import torch
import pandas as pd
import os


class PSSMDataset(torch.utils.data.Dataset):
    def __init__(self,pssm_in,label_in):
        self.pssm_lst = sorted(os.listdir(pssm_in))
        self.pssm_files = [pssm_in+f for f in self.pssm_lst]
        self.label_files = [label_in+f.split('.')[0]+'.label.txt' for f in self.pssm_lst]
    def __len__(self):
        return len(self.label_files)
    def __getitem__(self, idx):
        self.pssm_tensor = torch.tensor(pd.read_table(self.pssm_files[idx],index_col=0).astype(float).values)
        self.label_tensor = torch.tensor(list(map(int,open(self.label_files[idx],'r').readlines()[0].split(','))))
        return self.pssm_tensor, self.label_tensor

def my_collate(batch):
    data = [item[0] for item in batch]
    target = [item[1] for item in batch]
    dim_data = [d.size()[0] for d in data]
    max_seqlen = max(dim_data)
    data_transformed = [F.pad(d.T,(0,int(max_seqlen - d.size()[0])),'constant',0.0).T if d.size()[0]!=max_seqlen else d for d in data]
    target_transformed = [F.pad(d,(0,int(max_seqlen - d.size()[0])),'constant',0.0) if d.size()[0]!=max_seqlen else d for d in target]
    return [data_transformed, target_transformed]

pssminpath = '/Users/suhancho/data/Uniprot_metalbinding_challenge/processed_pssm/'
labelinpath = '/Users/suhancho/data/Uniprot_metalbinding_challenge/processed_label/'
train_ds = PSSMDataset(pssminpath,labelinpath)
train_dl = DataLoader(train_ds,batch_size = 4,shuffle = True,collate_fn=my_collate)

