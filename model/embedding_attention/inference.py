


from joblib import Parallel, delayed
inpath = '/Users/suhancho/data/Uniprot_metalbinding_challenge/neg_sequence/'
infiles = [inpath+f for f in os.listdir(inpath)]

Parallel(n_jobs=8)(delayed(inference_fasta)(i) for i in infiles)