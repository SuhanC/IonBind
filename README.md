# IonBind

## Materials and methods

Protein Sequence Data Processing

The protein sequences in POS_TRAIN_FULL.fasta were processed by chunking each sequence into segments of size 4, centered around each ion binding site. The same method was applied to the sequences in NEG_TRAIN.fasta, but with all positions in each sequence. The resulting windowed sequences were padded to a fixed length to enable batch processing. The chunked sequences containing metal ion binding sites were labeled with the corresponding metal ion, while sequences without metal binding ions were labeled as no binding metal ions. The labeled sequences were then transformed into embedding vectors and fed into the CNN model.

Model Structure and Training

The IonBind model was developed using Keras to predict ion binding sites using only protein sequence data (Figure 1). The model consisted of four convolutional layers and an attention layer, which learned the proximal context of the sequences around metal ion binding sites and non-binding sites. The fully connected layer mapped the prediction values to 30 classes of metal ions, including the non-binding class. The softmax function was applied to obtain the binding probability of each ion, and the class with the highest probability was compared to the real label to calculate the accuracy for each batch. The model performance was validated using 20% of the total dataset, and the model was trained for a total of 100 epochs with early stopping.

Tuning of Predicted Scores

The IonBind model was tuned using the prediction data. Before applying the argmax function to map the class of the sequence, the minimum value was set to test the performance improvement using the score thresholding method. The performance improvement was defined as the difference in the sum of evaluation metrics between the original prediction performance and the thresholded prediction performance. The minimum value for the argmax function was chosen as 0.6, which showed the highest improvement in performance (Figure 2,3).

Inference of Test Dataset

The final IonBind model was saved using model.save() and loaded to predict the test dataset. Althoughthe IonBind model was trained by batch, the test dataset contained a massive amount of raw sequencesthat could not be processed at once. Therefore, the inference was performed using parallel processing by generating multiple jobs using bash.




## Model structure used in the competition
![Model Structure](https://github.com/SuhanC/IonBind/blob/main/IonBind.png?raw=true)



### Reference papers

Sequence-based
> https://www.nature.com/articles/s41587-022-01307-0.pdf
> https://jcheminf.biomedcentral.com/articles/10.1186/s13321-022-00584-w
> https://www.biorxiv.org/content/10.1101/2022.04.27.489750v1.full
> https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6947422/
> https://www.nature.com/articles/s41598-021-03431-4#Sec14

Structure-based
> https://academic.oup.com/nar/article/50/W1/W13/6567476

Hybrid
> https://academic.oup.com/bioinformatics/article/36/10/3018/5753946
> https://github.com/Singh-Lab/dSPRINT
>

Representation
> https://ieeexplore.ieee.org/document/9477085
> https://github.com/sacdallago/bio_embeddings

Review
> https://reader.elsevier.com/reader/sd/pii/S2001037022002185?token=5E4D50B3BD3CA49D6CDE558EB6620782F302E16D521B5DC91FC469AB46C4D724EDE192A925B7F591D3577CF5B6D717E1&originRegion=us-east-1&originCreation=20220801025958

