def train_val_split(samples, split=0.2):
    RANDOM.shuffle(samples)
    n_val = int(len(samples) * split)
    return samples[:-n_val], samples[-n_val:]
    
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