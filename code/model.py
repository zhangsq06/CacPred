# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import torch
import torch.nn.functional as F
import torch.nn as nn
class BasicBlock(nn.Module):
	def __init__(self, in_planes, layers):
		super(BasicBlock, self).__init__()
		self.block = nn.Sequential(
			nn.BatchNorm2d(in_planes),
			nn.ReLU(),
			nn.Conv2d(in_planes, layers, (1,15), 1, (0,7)),
		)

	def forward(self, x):
		out = self.block(x)
		return torch.cat([x, out],1)
#######################################
class DenseBlock(nn.Module):
	def __init__(self, nb_layers, in_planes, grow_rate,):
		super(DenseBlock, self).__init__()
		layers = []
		for i in range(nb_layers):
			layers.append(BasicBlock(in_planes + i*grow_rate, grow_rate,))
		self.layer = nn.Sequential(*layers)
	def forward(self, x):
		return self.layer(x)     
######################################################
class model3(nn.Module):
    def __init__(self, ):
        super(model3, self).__init__()
        self.Conv1 = nn.Conv2d(in_channels=1, out_channels=32, kernel_size=[4,16])
        self.Conv2 = nn.Conv2d(in_channels=1, out_channels=32, kernel_size=[4,16])
        self.maxpool=nn.MaxPool2d(kernel_size=[1,16],stride=16)
        self.upsampl1 = nn.ConvTranspose2d(in_channels=32, out_channels=32, kernel_size=[1,16])
        self.upsampl2 = nn.ConvTranspose2d(in_channels=32, out_channels=32, kernel_size=[1,16])
        self.Linear1 = nn.Linear(88088,1)

        self.DenseBlock = DenseBlock(3,64,8)
        self.dropout=nn.Dropout(p=0.03)
    def forward(self, input1,input2):
        conv1=self.Conv1(input1)
        conv1 = F.relu(conv1)
        upsampl1=self.upsampl1(conv1)     
        #############################
        conv2 = self.Conv2(input2)
        conv2 = F.relu(conv2)
        upsampl2=self.upsampl2(conv2)
        
        cats = torch.cat((upsampl1,upsampl2),dim=1)
        dense1 = self.DenseBlock(cats)
        fc1=dense1.view(dense1.size(0),-1)
        fc1 = self.Linear1(fc1)
        fc1 = self.dropout(fc1)
        fc1=torch.squeeze(torch.sigmoid(fc1))
        return fc1     
