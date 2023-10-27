# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 19:23:14 2019

@author: shuangquanzhang
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
import util
import numpy as np
from compute_meme import obtain
import os
import argparse
#######################################
SEED = 1
torch.manual_seed(SEED)
torch.cuda.manual_seed(SEED)
##################inputs name####################
parser=argparse.ArgumentParser(description='demo')
parser.add_argument('--name',type=int,default='')
args=parser.parse_args()
####################
class MM(nn.Module):
    def __init__(self, ):
        super(MM, self).__init__()
        self.Conv1 = nn.Conv2d(1,32,(4,20))
    def forward(self, input):
        x = self.Conv1(input)
        x = F.relu(x) 
        return x
#######################       
def paddata(train_sequences):
    train_data_array = {}
    train_sequences_array = util.oneHot_data_encode(train_sequences)
    #train_sequences_array = util.add_appendix(train_sequences_array, 0.25, FILTER_LENGTH=9)
    alldata=train_sequences_array[:200]
    alldata=np.transpose(alldata,axes=[0,3,2,1])
    alldata=torch.from_numpy(alldata)
    return alldata

#########################
def obtain_motifs(ii):
    path='../model/'
    files=os.listdir(path)
    name = 'CEBPB_A549_CEBPB_Stanford_wgEncodeEH002818'
    train_sequences, train_sequences_alph,reverse_sequences,reverse_alph = util.load_data_encode(path_curr_data='../data/',data_path='../data/encode_1001/'+name[:]+'.seq.gz')
    print(train_sequences.shape,reverse_sequences.shape)
    traindata=paddata(train_sequences)
    reversedata=paddata(reverse_sequences)
    ##############################
    deep= MM()
    Model=torch.load(path+name+'.pkl')
    deep.weights=Model['Conv1.weight']
    deep.bias=Model['Conv1.bias']
    direct_signal=np.squeeze(deep(traindata).detach().numpy())
    ################################
    deep= MM()
    Model=torch.load(path+name+'.pkl')
    deep.weights=Model['Conv2.weight']
    deep.bias=Model['Conv2.bias']
    reverse_signal=np.squeeze(deep(reversedata).detach().numpy())

    print(direct_signal.shape,reverse_signal.shape)
    path='../motifs/'+name[:]
    if not os.path.exists(path):
        os.mkdir(path)
    #################################
    pp=0
    for j in range(direct_signal.shape[1]):
        seqname=path+'/filter_'+str(j)+'.txt'
        faname=path+'/filter_'+str(j)+'.fa'
        file1=open(seqname,'w')
        file2=open(faname,'w')
        for i in range(direct_signal.shape[0]):
            directindxx=direct_signal[i,j,:]
            directindex=directindxx.argmax()
            reverse_indxx=reverse_signal[i,j,:]
            reverseindex=reverse_indxx.argmax()
            if directindxx[directindex]>0  and directindex<=82 and reverse_indxx[i,j,reverseindex] >0:
                pp=1
                seq=str(train_sequences_alph[i]).strip('\n')[int(directindex):int(index+20)]
                file1.writelines(seq+'\n')
                strs='>'+'seq'+'_'+str(pp)+'_'+str(i+1)+'_'+str(directindex)+'_'+str(20)+'\n'
                file2.writelines(strs)
                file2.writelines(seq+'\n')              
        pp=0
######################################        
        file1.close()     
        file2.close()
    obtain(path,32)
if __name__=='__main__':
    obtain_motifs(ii=args.name)


