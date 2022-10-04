from preproc_utils import *
from sequence_postproc_utils import *
import pandas as pd
import numpy as np

def get_processed_pssm(pssm_file):
    test_pssm = get_PSSM_mat(pssm_file)
    test_seq = ''.join(test_pssm.index.tolist())
    test_pssm['Charge']=get_charge(test_seq)
    test_pssm['Hydro']=get_hydro(test_seq)

    test_pssm['Func_polar']=get_func_group(test_seq)[0]
    test_pssm['Func_alipathic']=get_func_group(test_seq)[1]
    test_pssm['Func_aromatic']=get_func_group(test_seq)[2]
    test_pssm['Func_negative']=get_func_group(test_seq)[3]
    test_pssm['Func_positive']=get_func_group(test_seq)[4]

    test_pssm['Volume_tiny']=get_volume(test_seq)[0]
    test_pssm['Volume_small']=get_volume(test_seq)[1]
    test_pssm['Volume_big']=get_volume(test_seq)[2]

    test_pssm['HBdonor']=get_h_bond(test_seq)[0]
    test_pssm['HBacceptor']=get_h_bond(test_seq)[1]
    return(test_pssm)
 


 