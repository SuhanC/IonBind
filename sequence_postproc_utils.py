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
    