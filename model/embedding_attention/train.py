from csv import DictReader
from random import Random
from collections import Counter
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from tqdm import tqdm
from keras.utils.np_utils import to_categorical
import pickle as pkl
import warnings
warnings.filterwarnings('ignore')
import os
from sklearn.metrics import classification_report


RANDOM = Random(42)

dat_prot = pd.read_table('/Users/suhancho/data/Uniprot_metalbinding_challenge/sequence_df.tsv')
chebi = pd.read_table('/Users/suhancho/data/Uniprot_metalbinding_challenge/POS_TRAIN_FULL.tsv')
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

vocabulary = build_vocabulary(train_samples)

train_samples = [transform(sample, vocabulary) for sample in train_samples]
val_samples = [transform(sample, vocabulary) for sample in val_samples]



model = predict_bindsite(
    vocab_size=len(vocabulary),
    char_embedding_size=16,
    base_filters=32,#32
    doc_embedding_size=300,
    dropout=0.1
)
model.summary()



train_batch_provider = BatchProvider(train_samples, batch_size=4096, shuffle=True, run_forever=True)
train_batches = train_batch_provider.generate_batches()
val_batch_provider = BatchProvider(val_samples, batch_size=4096, shuffle=False, run_forever=True)
val_batches = val_batch_provider.generate_batches()


from keras.callbacks import ModelCheckpoint, EarlyStopping
EPOCHS = 100
model_checkpoint_callback = ModelCheckpoint('./weights.hdf5',monitor='val_acc',save_best_only=True)
early_stopping_callback = EarlyStopping(monitor='val_loss', patience=5)

model.fit_generator(
        generator=train_batches,
        steps_per_epoch=len(train_batch_provider),
        epochs=EPOCHS,
        validation_data=val_batches,
        validation_steps=len(val_batch_provider),
        callbacks= [model_checkpoint_callback,early_stopping_callback]
    )


model.save('/Users/suhancho/script/uniprot/IonBind/saved_models/Ionbind_221214_9mer_epoch100')
model = keras.models.load_model("/Users/suhancho/script/uniprot/IonBind/saved_models/Ionbind_221214_9mer_epoch100")




val_batch_provider = BatchProvider(val_samples, batch_size=2056, shuffle=False, run_forever=False)
val_predictions = model.predict_generator(val_batch_provider.generate_batches(), steps=len(val_batch_provider))
val_targets = np.array([sample['target'] for sample in val_samples])
val_prediction_argmax = [np.argmax(p) for p in val_predictions]
val_targets_argmax = [np.argmax(p) for p in val_targets]

report = classification_report(val_targets_argmax,val_prediction_argmax,output_dict = True)
report = pd.DataFrame(report).T
# Target ion information processing and result processing
target_ion_info = dat_prot_bindsite_sampled[['target','ChEBI-ID']]
target_ion_info['target'] = target_ion_info['target'].apply(np.argmax).astype(int)
target_ion_info.drop_duplicates(inplace = True)
report = classification_report(val_targets_argmax,val_prediction_argmax,output_dict = True)
report = pd.DataFrame(report)
report_class = report.drop(['accuracy','macro avg','weighted avg'],axis=1).T
report_class = report_class.reset_index()
report_class['index'] = report_class['index'].astype(int)
report_class = pd.merge(report_class,target_ion_info,left_on = 'index',right_on = 'target').drop('index',axis=1)
report_class = report_class[['target','ChEBI-ID','precision','recall','f1-score','support']].T

report_class.to_csv('./Submission.csv',index = None)