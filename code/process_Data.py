# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 21:05:21 2019

@author: 11154
"""
import re
import os
import numpy as np
import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--name',default='TEST_peaks.bed',help='bed')
parser.add_argument('--flank',default=500,type=int)
parser.add_argument('--genome',default='hg19',type=str)
args = parser.parse_args()
################################################################
def get_peaks(args):
    name=args.name
    name=name.split('.')[0]
    file=open(args.name,'r')
    data=file.readlines() 
    file.close()
    name='../TfbsUniform_hg19_ENCODE/'+name+'.narrowPeak'
    file1=open(name,'w')
    pattern=re.compile(r'chr')
    for t in data:
        match_res = pattern.match(t)
        if match_res:
            t1=t.split('\t')
            tt=str(t1[0])+'\t'+ str(t1[1]) + '\t' + str(t1[2])+'\n'
            file1.writelines(tt) 
    file1.close()
    cmd='gzip '+ name
    os.system(cmd)
#############################################################
def getchrom_index(args):
    name=args.name
    name=name.split('.')[0]
    file=open(args.name,'r')
    data1=file.readlines() 
    file.close()
    sortn=[]
    ll=len(data1)
    for i in range(ll):
        t=data1[i]
        t1=t.strip('\n').split('\t')[4]
        sortn.append(t1)
    sortn=np.array(sortn,dtype=np.float32)
    sortn_index=sortn.argsort()[::-1]
    test_val=[name+'_train.bed',name+'_test.bed']
    test_val_fa=[name+'train1.fa',name+'test1.fa']
    pattern=re.compile(r'chr')
    file0=open(test_val[0],'w')
    file1=open(test_val[1],'w')
    for j in range(ll):
        t=data1[sortn_index[j]]
        match_res = pattern.match(t)
        t1=t.split('\t')
        if match_res:
            t11=round((int(t1[1])+int(t1[2]))/2)-int(args.flank)+1
            t12=round((int(t1[1])+int(t1[2]))/2)+int(args.flank)+2
            if t11 >0 and t12 > 0:
                tt=str(t1[0])+'\t'+ str(int(t11)) + '\t' + str(int(t12))+'\t'+str(t1[4])
                if j % 2==0:
                    file0.writelines(tt) 
                else:
                    file1.writelines(tt)      
    file0.close()
    file1.close()
    for k in range(2):
        cmd='bedtools getfasta -fi ../GRCh37.p13.genome.fa -bed '+ test_val[k]+ ' -s -fo '+test_val_fa[k]
        os.system(cmd)
def readbeds(name):
    file=open(name,'r')
    data=file.readlines()
    file.close()
    return data
def makeseq(input_namefa,input_namescore,out_name):
    data=readbeds(input_namefa)
    datascores=readbeds(input_namescore)
    data1=open(out_name,'w')
    times=0
    for t1 in data:
        if t1.startswith('>'):
            for score in datascores:
                starts=score.split('\t')[1]
                stops=score.split('\t')[2]
                if (starts in t1) and (stops in t1):
                    scs=score.split('\t')[3]       
        else:
            ttt='A'+'\t'+'peaks'+'\t'+t1.strip().upper()+'\t'+str(scs)
            data1.writelines(ttt) 
            times+=1
    data1.close() 
######################################################################
def get_data(args):
    name=args.name
    name=name.split('.')[0]
    innames=[name+'train1.fa',name+'test1.fa']
    test_val=[name+'_train.bed',name+'_test.bed']
    outnames=[name+'_AC.seq',name+'_B.seq',name+'.seq']
    for i in range(2):
        makeseq(innames[i],test_val[i],outnames[i])
        cmd='cat '+outnames[i]+' >> '+outnames[2]
        os.system(cmd)
        cmd='gzip '+ outnames[i]
        os.system(cmd)
        cmd='mv '+ outnames[i]+'.gz'+' ../encode_'+str(args.flank*2+1)
        os.system(cmd)
#        os.remove(innames[i])
    cmd='gzip '+ outnames[2]
    os.system(cmd)
    cmd='mv '+ outnames[2]+'.gz'+' ../encode_'+str(args.flank*2+1)
    os.system(cmd)
    file=open('../encode_tfbs.txt','a+')
    txxt=name+'_encode'+'\t'+name+''+'\n'
    file.writelines(txxt)
    file.close()
    files1=open(innames[0],'r')
    datas1=files1.readlines()
    files1.close()
    files2=open(innames[1],'r')
    datas2=files2.readlines()
    files2.close()
    if len(datas1)>len(datas2):
        ll=len(datas2)
    else:
        ll=len(datas1)
    allfile=open(name+'_all.fa','w')
    for i in range(0,ll-2,2):
        allfile.writelines(datas1[i])
        allfile.writelines(datas1[i+1])
        allfile.writelines(datas2[i])
        allfile.writelines(datas2[i+1])
    allfile.close()
    os.remove(innames[0])
    os.remove(innames[1])
#################get_bed######################################
    file=open(args.name,'r')
    data1=file.readlines()
    file.close()
    out=open(name+'_all.bed','w')
    ll=len(data1)
    for i in range(ll):
        t=data1[i]
        t1=t.strip('\n').split('\t')
        t11=round((int(t1[1])+int(t1[2]))/2)-int(args.flank)+1
        t12=round((int(t1[1])+int(t1[2]))/2)+int(args.flank)+2
        tt=str(t1[0])+'\t'+ str(int(t11)) + '\t' + str(int(t12)) +'\n'
        out.writelines(tt)
    out.close()
if __name__=='__main__':
    get_peaks(args)
    getchrom_index(args)
    get_data(args)
                
