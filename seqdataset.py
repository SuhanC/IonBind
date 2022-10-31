import pandas as pd
import os
import numpy as np


def get_charge(sequence):
    transformed_seq=[]
    positive=['K','R']
    negative=['D','E']    
    for s in sequence:
        if positive.count(s):
            transformed_seq+=[1]
        elif negative.count(s):
            transformed_seq+=[-1]
        else : 
            transformed_seq+=[0]
    return(transformed_seq)

def get_func_group(sequence):
    polar=['C','N','P','Q','S','T','U']
    alipathic=['A','G','I','L','M','V','J']
    aromatic=['F','W','Y']
    negative=['D','E']
    positive=['H','K','R']

    def get_func_representation(sequence,group_idx):
        return([1 if group_idx.count(s) else 0 for s in sequence])
        
    result_lst=[]
    for func_group in [polar,alipathic,aromatic,negative,positive]:
        result_lst.append(get_func_representation(sequence,func_group))
    return(result_lst)


def get_hydro(sequence):
    transformed_seq=[]
    hydrophilic=['D','E','G','K','N','P','Q','R','S','T','W','Y','Z']
    for s in sequence:
        if hydrophilic.count(s):
            transformed_seq+=[1]
        else : 
            transformed_seq+=[0]
    return(transformed_seq)

def get_volume(sequence):
    tiny = ['A','C','G','S','U']
    small = ['D','N','P','T','V','B']
    big = ['E','F','H','I','K','L','M','Q','R','W','Y','Z','J']
    def get_func_representation(sequence,group_idx):
        return([1 if group_idx.count(s) else 0 for s in sequence])
    result_lst=[]
    for volume_group in [tiny,small,big]:
        result_lst.append(get_func_representation(sequence,volume_group))
    return(result_lst)


def get_h_bond(sequence):   
    aa_h_bond_donor = {
        'A': 0.95,
        'C': 0.96,
        'D': 0.89,
        'E': 0.93,
        'F': 0.92,
        'G': 0.96,
        'H': 2.62,
        'I': 0.93,
        'K': 3.59,
        'L': 0.98,
        'M': 1.02,
        'N': 2.24,
        'P': 0.0,
        'Q': 2.27,
        'R': 5.4,
        'S': 2.31,
        'T': 2.29,
        'V': 0.93,
        'W': 1.83,
        'Y': 2.00,
        'X': 0.0,
        '*': 0.0,
        'U': 0.0,
        'B': 1.5,
        'J': 0.95,
        'Z': 1.5

    }

    aa_h_bond_acceptor = {
        'A': 1.16,
        'C': 1.12,
        'D': 4.84,
        'E': 4.43,
        'F': 1.08,
        'G': 1.18,
        'H': 1.94,
        'I': 1.12,
        'K': 1.12,
        'L': 1.20,
        'M': 1.12,
        'N': 2.41,
        'P': 1.13,
        'Q': 2.41,
        'R': 1.20,
        'S': 2.04,
        'T': 2.01,
        'V': 1.11,
        'W': 1.21,
        'Y': 1.68,
        'X': 0.0,
        '*': 0.0,
        'U': 0.0,
        'B': 3.50,
        'J': 1.15,
        'Z': 3.00


    }
    result_donor=[aa_h_bond_donor[s] for s in sequence]
    result_acceptor=[aa_h_bond_acceptor[s] for s in sequence]
    return([result_donor,result_acceptor])
    
def get_PSSM_mat(pssm_file):
    test = open(pssm_file,'r').readlines()

    aa_cols=['Residue','A', ' R', ' N', ' D', ' C', ' Q', ' E', ' G', ' H', ' I', ' L', ' K', ' M', ' F', ' P', ' S', ' T', ' W', ' Y', ' V']

    pssm_lst=[]
    for i in range(3, len(test)):
        pssm_line = test[i].split('  ')
        pssm_line = [l.strip() for l in pssm_line if l!='']
        matrix_1 = pssm_line[0:21]
        matrix_1 = [matrix_1[i].split(' ')[0]if i!=0 else matrix_1[i] for i in range(len(matrix_1))]
        matrix_2 = pssm_line[21:41]
        pssm_lst.append(matrix_1)

    pssm = pd.DataFrame(pssm_lst)
    pssm.columns = aa_cols
    pssm = pssm.dropna()

    pssm['Residue'] =[r.split(' ')[1] for r in pssm['Residue'].tolist()]
    pssm.set_index('Residue',inplace = True)
    return(pssm)


def get_processed_pssm(pssm_mat):
    test_seq = ''.join(pssm_mat.index.tolist())
    pssm_mat['Charge']=get_charge(test_seq)
    pssm_mat['Hydro']=get_hydro(test_seq)

    pssm_mat['Func_polar']=get_func_group(test_seq)[0]
    pssm_mat['Func_alipathic']=get_func_group(test_seq)[1]
    pssm_mat['Func_aromatic']=get_func_group(test_seq)[2]
    pssm_mat['Func_negative']=get_func_group(test_seq)[3]
    pssm_mat['Func_positive']=get_func_group(test_seq)[4]

    pssm_mat['Volume_tiny']=get_volume(test_seq)[0]
    pssm_mat['Volume_small']=get_volume(test_seq)[1]
    pssm_mat['Volume_big']=get_volume(test_seq)[2]

    pssm_mat['HBdonor']=get_h_bond(test_seq)[0]
    pssm_mat['HBacceptor']=get_h_bond(test_seq)[1]
    return(pssm_mat)

def process_binddata(bind_filepath):
    bind_tsv_list = [bind_filepath + f for f in os.listdir(bind_filepath)]
    bindlist = pd.concat([pd.read_table(f) for f in bind_tsv_list])
    bindlist['Name'] = pd.Categorical(bindlist.Name)
    bindlist['Target'] = bindlist.Name.cat.codes
    return(bindlist)

def get_label(accession,infile,mode,bind_tsv):
    if mode=='positive':
        sequence_fasta = open(infile,'r').readlines()[0]
        bindset = bind_tsv[bind_tsv.Accession==accession].sort_values('Position')

        if len(bindset)!=0:
            label=['0']*len(sequence_fasta)
            for i in range(len(bindset)):
                pos = bindset.Position.tolist()[i]-1
                label[pos] = bindset.Target.tolist()[i]
        else : 
            label = np.NaN
    elif mode=='negative':
        sequence_fasta = open(infile,'r').readlines()[0]
        label = ['29']*len(sequence_fasta)
    return(label)