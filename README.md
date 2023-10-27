# CacPred: a cascaded convolutional neural network for TF-DNA binding prediction 
# Dependencies
Python3.6, Biopython-1.72, pytorch-1.9.1, meme 4.5.0, bedtools v2.21.0. 
# Data
Download GRCh37.p13.genome.fa from https://www.gencodegenes.org/human/release_19.html, and put it into the current directory.
The first step make the 'code', ‘TfbsUniform_hg19_ENCODE’, ‘encode_1001’, ’model’, 'motifs' and ‘output’ directory. And then move all python scripts into the 'code' directory.
# A simple tutorial
## 1 Preprocess the bed file.

__Usage:__ python process_Data.py --bed  TEST_peaks.bed <br>
Arguments: 
--bed (the prefix name of your file).  

## 2 Training and testing a CacPred model 
__Usage:__ python train_test.py --index 0  <br>
Arguments:   
-- index (the index of the dataset in encode_tfbs.txt).  
When finishing the process of training and testing, the predicted results are stored in the ‘output’ dir, and trained models are stored in the ‘model’ dir. 
## 3 Finding motifs
__Usage:__ python detect_motifs.py –-name TEST_peaks <br>
Arguments:  
--bed (the prefix name of your file).  
This command will output the found motifs to the ‘motifs’ dir.
