# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 12:02:22 2018

@author: 11154
"""

from sklearn import metrics  
def com_auc(test_l,pred):
    ## 真实值，预测值
    fpr, tpr, thresholds = metrics.roc_curve(test_l, pred, pos_label=1)
    roc_auc = metrics.auc(fpr, tpr)
    return roc_auc,fpr,tpr,thresholds
 
  
    