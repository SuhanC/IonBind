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

    # inference_result.append([testprot_og,test_prediction_argmax,test_prediction_proba])
    # return([testprot_og,test_prediction_argmax,test_prediction_proba])
    with open(outpath+testprot_og+'.result.txt','w') as outfile:
        outfile.write('Protein_Name\tIon\tPrediction_Score\n')
        outfile.write('\t'.join([testprot_og,','.join(test_prediction_argmax),','.join(test_prediction_proba)]))
        outfile.close()


        
from joblib import Parallel, delayed
inpath = '/Users/suhancho/data/Uniprot_metalbinding_challenge/neg_sequence/'
infiles = [inpath+f for f in os.listdir(inpath)]

Parallel(n_jobs=8)(delayed(inference_fasta)(i) for i in infiles)