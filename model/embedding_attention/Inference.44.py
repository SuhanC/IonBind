from csv import DictReader
from random import Random
from collections import Counter
import numpy as np
import pandas as pd
from tqdm import tqdm
from keras.utils.np_utils import to_categorical
import pickle as pkl
import warnings
warnings.filterwarnings('ignore')
import os
import keras
from keras.utils import GeneratorEnqueuer
import csv
import gc
import tensorflow as tf

RANDOM = Random(42)

def load_samples(path):
    with open(path) as f:
        samples = list(DictReader(f))
        for sample in samples:
            sample['target'] = int(sample.get('target', -1))
        return samples


def train_val_split(samples, split=0.2):
    RANDOM.shuffle(samples)
    n_val = int(len(samples) * split)
    return samples[:-n_val], samples[-n_val:]


import pandas as pd
dat_prot = pd.read_table('/home/suhan/data/Uniprot_metalbinding_challenge/sequence_df.tsv')
chebi = pd.read_table('/home/suhan/data/Uniprot_metalbinding_challenge/POS_TRAIN_FULL.tsv')

dat_prot_bindsite = pd.merge(dat_prot,chebi,left_on = 'Protein name',right_on='Accession',how = 'outer')
neg_seqs = dat_prot_bindsite[dat_prot_bindsite.Type=='neg_sequence']
neg_seqs['Position'] = [int(divmod(len(p),2)[0]) for p in neg_seqs['Protein sequence'].tolist()]
neg_seqs['ChEBI-ID'] = 'NB'
dat_prot_bindsite = pd.concat([dat_prot_bindsite[dat_prot_bindsite['Type']=='pos_sequence'],neg_seqs]).reset_index(drop = True)
dat_prot_bindsite['ChEBI-ID'] = pd.Categorical(dat_prot_bindsite['ChEBI-ID'])
dat_prot_bindsite['target'] = dat_prot_bindsite['ChEBI-ID'].cat.codes
dat_prot_bindsite['Position'] = dat_prot_bindsite['Position'].astype(int)
dat_prot_bindsite_sampled = dat_prot_bindsite.sample(frac=1,random_state=9510)
# dat_prot_bindsite_sampled.to_csv('/Users/suhancho/data/Uniprot_metalbinding_challenge/data_before_windowing.tsv',sep='\t')


def get_windowdf(protfile,winsize=4):
    window=[]
    with open(protfile,'r') as p:
        for line in tqdm(p):
            if not line.count('Protein sequence'):
                protseq = line.split('\t')[2].strip()
                bindsite = int(line.split('\t')[-2].strip())

                if ((bindsite-winsize)<0 )& ((bindsite+winsize) > len(protseq)) :
                    front_pad_seq = 'X'*(winsize-bindsite)
                    bindseq = protseq[bindsite-winsize : bindsite+winsize]
                    rear_pad_seq = 'X'*((winsize+bindsite)-len(protseq))
                    windowed = front_pad_seq+bindseq+rear_pad_seq

                elif ((bindsite-winsize)<0) & ((bindsite+winsize) < len(protseq)):
                    front_pad_seq = 'X'*(winsize-bindsite)
                    bindseq = protseq[bindsite-winsize : bindsite+winsize]
                    windowed = front_pad_seq+bindseq

                elif ((bindsite-winsize)>0) & ((bindsite+winsize) > len(protseq)):
                    rear_pad_seq = 'X'*((winsize+bindsite)-len(protseq))
                    bindseq = protseq[bindsite-winsize : bindsite+winsize]
                    windowed = bindseq+rear_pad_seq

                else:
                    windowed = protseq[bindsite-winsize : bindsite+winsize]
                    
                window.append(windowed)
    return(window)

def build_vocabulary(samples, vocab_min_freq=100):
    counts = Counter(ch for sample in samples for ch in sample['question_text'])
    chars = sorted(ch for ch, count in counts.items() if count >= vocab_min_freq)
    return {char: i for i, char in enumerate(chars)}

def transform(sample, vocabulary):
    sample['encoded_text'] = [vocabulary[ch] for ch in sample['question_text'] if ch in vocabulary]
    return sample


test = get_windowdf('/home/suhan/data/Uniprot_metalbinding_challenge/data_before_windowing.tsv',winsize = 4)
dat_prot_bindsite_sampled['window_4'] = test
dat_prot_bindsite_sampled = dat_prot_bindsite_sampled.reset_index()
dat_prot_bindsite_sampled.rename(columns = {'index':'qid','window_4':'question_text'},inplace = True)
dat_prot_bindsite_sampled['qid'] = 'HASH_'+dat_prot_bindsite_sampled['qid'].astype(str)
one_hot_labels = to_categorical(dat_prot_bindsite_sampled['target'].tolist())
dat_prot_bindsite_sampled['target'] = list(one_hot_labels)
samples=dat_prot_bindsite_sampled[['qid','question_text','target']].to_dict('records')
train_samples, val_samples = train_val_split(samples)

vocabulary = build_vocabulary(train_samples)

train_samples = [transform(sample, vocabulary) for sample in train_samples]
val_samples = [transform(sample, vocabulary) for sample in val_samples]

from math import ceil

class TestBatchProvider:
    def __init__(self, samples, batch_size, shuffle=False, run_forever=False):
        self._samples = samples
        self._batch_size = batch_size
        self._shuffle = shuffle
        self._run_forever = run_forever

    def generate_batches(self):
        batch = []
        indices = list(range(len(self._samples)))
        while True:
            if self._shuffle:
                RANDOM.shuffle(indices)
            for i in indices:
                batch.append(self.get_item(i))
                if len(batch) == self._batch_size:
                    yield self.transform_batch(batch)
                    batch = []
            if not self._run_forever:
                break
        if batch:
            yield self.transform_batch(batch)

    def __len__(self):
        return int(ceil(len(self._samples) / self._batch_size))

    def get_item(self, idx):
        sample = self._samples[idx]
        return sample['encoded_text']

    def transform_batch(self, items):
        texts= items
        max_length = max(len(text) for text in texts)
        text_batch = np.zeros((len(texts), max_length))
        for i, text in enumerate(texts):
            text_batch[i, :len(text)] = text
        return text_batch


def inference_fasta(infile):
    model = keras.models.load_model("Saved_Ionbind_NLP_221210")

    outpath='/home/suhan/data/Uniprot_metalbinding_challenge/inference_result/'
    testprot_tmp = infile.split('/')[-1].split('.')[0]
    testprot_og = infile.split('/')[-1].split('.')[0]
    testseq_tmp_file = open(infile,'r')
    testseq_tmp = testseq_tmp_file.readlines()[1].strip()
    test_input_lst=[]
    for i in list(range(len(testseq_tmp))):
        test_input={} # Generate Dictionary for prediciton 
        testprot_tmp = 'HASH_'+str(testprot_og)+'.'+str(i)
        if i-4<0:
            windowed_tmp = 'X'*(4-i)+testseq_tmp[0:i+5]
        elif i+5>len(testprot_tmp):
            windowed_tmp = testseq_tmp[i-4:i+5]+'X'*(i+5-len(testseq_tmp)+1)
        else : 
            windowed_tmp = testseq_tmp[i-4:i+5]

        test_input['qid'] = testprot_tmp
        test_input['question_text'] = windowed_tmp
        test_input_lst.append(test_input)

    test_input_lst = [transform(sample, vocabulary) for sample in test_input_lst]
    test_batch_provider = TestBatchProvider(test_input_lst, batch_size=len(test_input_lst), shuffle=False, run_forever=False)
    enqueuer = GeneratorEnqueuer(test_batch_provider.generate_batches())
    enqueuer.start()
    test_batches = enqueuer.get()
    for batch in test_batches:
        test_predictions = model.predict_on_batch(batch)
        test_prediction_argmax = [np.argmax(p) for p in test_predictions]
        test_prediction_proba = [str(round(prob[idx],4)) for prob,idx in zip(test_predictions,test_prediction_argmax)]
        test_prediction_argmax = [str(p) for p in test_prediction_argmax]
        tf.keras.backend.clear_session()
        gc.collect()
        testseq_tmp_file.close()
    with open(outpath+testprot_og+'.result.txt','w') as outfile:
        outfile.write('Protein_Name\tIon\tPrediction_Score\n')
        outfile.write('\t'.join([testprot_og,','.join(test_prediction_argmax),','.join(test_prediction_proba)]))
        outfile.close()
    tf.keras.backend.clear_session()
    gc.collect()
    return(0)

#################################Run###############################

from joblib import Parallel, delayed
# inpath = '/home/suhan/data/Uniprot_metalbinding_challenge/test_sequence/'
# infiles = [inpath+f for f in os.listdir(inpath)]
infiles = open('/home/suhan/data/Uniprot_metalbinding_challenge/testfiles.txt','r').readlines()
infiles = [f.strip('\n') for f in infiles]

Parallel(n_jobs=64, backend='loky')(delayed(inference_fasta)(i) for i in infiles)