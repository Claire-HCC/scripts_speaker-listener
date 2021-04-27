# coding: utf-8
# #!/bin/env python
import brainiak.eventseg.event
import numpy as np
import glob
import os
import scipy.io
from scipy.stats import stats
import logging
import matplotlib.pyplot as plt
from scipy.stats import norm, zscore, pearsonr, ttest_1samp
import pickle
import time
import sys
import re

# this should finish in 3hr for one roi

print(sys.argv)
exp=(sys.argv[1])
print(exp)
ri=int(sys.argv[2])
print(ri)

expdir = '/scratch/claire/speaker-listener/'
timeUnit='tr'
# load eventN_test
eventN_test=range(10,121)
w=5
froidir='mor'
                             
start_time = time.time()

# load mat file
fnames=glob.glob(os.path.join(expdir + '/' + exp + '/fmri/timeseries/' + timeUnit + '/network/' + froidir , 'zscore_listenerAll_iscmasked_*.mat'))
rname=re.search('/scratch/claire/speaker-listener//'+exp+'/fmri/timeseries/tr/network/mor/zscore_listenerAll_iscmasked_(.+?).mat', fnames[ri-1]).group(1) 
data_mat=scipy.io.loadmat(expdir+ exp + "/fmri/timeseries/tr/network/"+froidir+"/zscore_listenerAll_" + rname +".mat")
gdata=data_mat['gdata']

tn=gdata.shape[1]
subjn=gdata.shape[2]
voxn=gdata.shape[0]

within_across_real = np.empty([np.max(eventN_test),subjn])    
within_across_real[:]=np.nan

for K in eventN_test:
    Ki=K-1
    
    d=np.load(expdir +  exp+ '/fmri/hmm/ListenersLeave1Out_eventLabels/'+ rname + '_ListenersLeave1Out_eventLabels_K' +str(K)+ '.npz')
    eventLabels_self=d['eventLabels_self']
    
    for s in range(subjn):
        if np.sum(~np.isnan(gdata[:,:,s]))>0:
            self=gdata[:,:,s]
        
        corrs = np.zeros(tn-w)
        for t in range(tn-w):
            corrs[t] = pearsonr(self[:,t],self[:,t+w])[0]
        
        eventLabels=eventLabels_self[s,:]
        _, event_lengths = np.unique(eventLabels, return_counts=True)
        
        within = corrs[eventLabels[:-w] == eventLabels[w:]].mean()
        across = corrs[eventLabels[:-w] != eventLabels[w:]].mean()
        within_across_real[Ki,s] = within - across

np.save(expdir +  exp+ '/fmri/hmm/WithinAcrossSim/'+ rname + '_ListenersLeave1Out_withinAcross.npy', within_across_real)      
print(time.time()-start_time)    
