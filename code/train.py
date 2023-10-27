# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 15:01:55 2019

@author: shuangquanzhang
"""
import torch
import torch.optim as optim
import numpy as np
import torch.nn as nn
import torch.utils.data as Data
from util import convert_data
from auc import com_auc
from model import model3
import argparse
torch.manual_seed(1)
np.random.seed(1)
torch.cuda.current_device()
torch.cuda._initialized = True
parser = argparse.ArgumentParser(description='input some paras.')
parser.add_argument('--index',default=0,type=int)
parser.add_argument('--m',default='model3',type=str)
args = parser.parse_args()
EPOCH = 20
LR = 5e-2
##########################
traindirectx, traindirecty, trainreversex, name = convert_data(args.index,flag='train')
traindirectx, trainreversex = traindirectx[:,np.newaxis,:,:], trainreversex[:,np.newaxis,:,:]
directy=np.squeeze(directy)
nums = len(directy)
validX1_data = torch.FloatTensor(traindirectx[int(nums*0.8):])
validX2_data = torch.FloatTensor(trainreversex[int(nums*0.8):])        
validY_data = torch.FloatTensor(traindirecty[int(nums*0.8):])

trainX1_data = torch.FloatTensor(traindirectx[:int(nums*0.8)])
trainX2_data = torch.FloatTensor(trainreversex[:int(nums*0.8)])
trainY_data = torch.FloatTensor(traindirecty[int(nums*0.8):])
########################################
params = {'batch_size': 100}
train_loader = Data.DataLoader(
    dataset=Data.TensorDataset(trainX1_data, trainX2_data, trainY_data), 
    shuffle=False,
    **params)
print('compling the network')
##########################################
if args.m=='model2':
    MODEL=model2()
if args.m=='model3':
    MODEL=model3()
###################training data##########
MODEL.cuda()
optimizer = optim.Adadelta(MODEL.parameters(), lr=LR)
scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=3,verbose=True)
loss_func = nn.BCELoss()        
print('starting training')
val_loss=10        
for epoch in range(EPOCH):
    MODEL.train(mode=True)
    for step, (direct_batch_x, reverse_batch_x, train_batch_y) in enumerate(train_loader):       
        direct_batch_x = direct_batch_x.cuda()
        reverse_batch_x = reverse_batch_x.cuda()
        train_batch_y = train_batch_y.cuda()
        out = MODEL(direct_batch_x, reverse_batch_x)
        loss = loss_func(out, train_batch_y)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
#########################################
    if epoch % 10 ==0:
        with torch.no_grad():
            for i,j in zip(range(0,validX1_data.shape[0],100), range(100,validX2_data.shape[0]+100,100)):
                valout = MODEL(validX1_data[i:j].cuda(), validX2_data[i:j].cuda())
                valout=valout
                if i==0:
                    val_out=valout
                else:
                    val_out=torch.cat((val_out,valout),0)      
            valid_loss = loss_func(val_out, validY_data.cuda())
            scheduler.step(valid_loss)
            val_auc, _, _, _ = com_auc(validY_data.data.cpu().numpy(),val_out.data.cpu().numpy())
            print('valid_loss: %f ' % valid_loss.data.cpu().numpy())
#################################
model_name='../model/'+name+'.pkl'
torch.save(MODEL.state_dict(), model_name) 
##################Testing model################
print('.....testing.........')
testX1_data, testX2_data, testY_data, _ = convert_data(args.index, flag='test')
testX1_data, testX2_data = testX1_data[:,np.newaxis,:,:], testX2_data[:,np.newaxis,:,:]
testY_data = np.squeeze(testY_data)
testX1_data = torch.FloatTensor(testX1_data)
testX2_data = torch.FloatTensor(testX2_data)
##################################################################################################
with torch.no_grad():
    for i,j in zip(range(0, testX1_data.shape[0],100), range(100, testX2_data.shape[0]+100,100)):
        testout = MODEL(testX1_data[i:j].cuda(), testX2_data[i:j].cuda())
        testout=testout.data.cpu().numpy()
        if i==0:
            test_out=testout
        else:
            test_out=np.concatenate((test_out,testout),0)
        testout=[]
################################################################################################## 
test_auc,_,_,_ = com_auc(np.array(testY_data), test_out)
print('test_auc: %f ' % test_auc)
out_name='../output/'+name+'.txt'
np.savetxt(out_name,[np.array(testY_data), test_out])
