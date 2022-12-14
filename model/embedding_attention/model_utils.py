from math import ceil
import gc
import tensorflow as tf
from keras.utils import GeneratorEnqueuer

class BatchProvider:
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
        return sample['encoded_text'], sample['target']

    def transform_batch(self, items):
        texts, targets = zip(*items)
        max_length = max(len(text) for text in texts)
        text_batch = np.zeros((len(texts), max_length))
        for i, text in enumerate(texts):
            text_batch[i, :len(text)] = text
        target_batch = np.array(targets)
        return text_batch, target_batch


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
    outpath='/Users/suhancho/data/Uniprot_metalbinding_challenge/inference_result/'
    testprot_tmp = infile.split('/')[-1].split('.')[0]
    testprot_og = infile.split('/')[-1].split('.')[0]
    testseq_tmp_file = open(infile,'r')
    testseq_tmp = testseq_tmp_file.readlines()[0].strip()
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
        outfile.write('|'.join([testprot_og,','.join(test_prediction_argmax),','.join(test_prediction_proba)]))
        outfile.close()
