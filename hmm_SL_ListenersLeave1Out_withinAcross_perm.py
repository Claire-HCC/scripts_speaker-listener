# coding: utf-8
# #!/bin/env python
import brainiak.eventseg.event
import numpy as np
import glob
import os
import scipy.io
import logging
import time
import sys
import re
from scipy.stats import norm, zscore, pearsonr, ttest_1samp
# this should finish in 3hr for one roi

print(sys.argv)
exp=(sys.argv[1])
print(exp)
ri=int(sys.argv[2])
print(ri)

expdir = '/scratch/claire/speaker-listener/'
timeUnit='tr'
exps=['pieman','bronx','merlin','sherlock']
# load eventN_test
eventN_test=range(10,121)
w=5
nPerm=1000
froidir='mor'

start_time = time.time()

# load mat file
fnames=glob.glob(os.path.join(expdir + '/' + exp + '/fmri/timeseries/' + timeUnit + '/roi/' + froidir , 'zscore_listenerAll_*.mat'))
rname=re.search('/scratch/claire/speaker-listener//merlin/fmri/timeseries/tr/roi/mor/zscore_listenerAll_(.+?).mat', fnames[ri-1]).group(1) 
data_mat=scipy.io.loadmat(expdir+ exp + "/fmri/timeseries/tr/roi/"+froidir+"/zscore_listenerAll_" + rname +".mat")
gdata=data_mat['gdata']

tn=gdata.shape[1]
subjn=gdata.shape[2]
voxn=gdata.shape[0]

within_across_null = np.empty([np.max(eventN_test),subjn, nPerm])
within_across_null[:]=np.nan

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
        
        np.random.seed(0)
        for p in range(nPerm+1):
            perm_lengths = np.random.permutation(event_lengths)
            eventLabels = np.zeros(tn, dtype=np.int)
            eventLabels[np.cumsum(perm_lengths[:-1])] = 1
            eventLabels = np.cumsum(eventLabels)
            
            within = corrs[eventLabels[:-w] == eventLabels[w:]].mean()
            across = corrs[eventLabels[:-w] != eventLabels[w:]].mean()
            within_across_null[Ki,s,p] = within - across

np.save(expdir +  exp+ '/fmri/hmm/WithinAcrossSim/'+ rname + '_ListenersLeave1Out_withinAcross_perm.npy', within_across_null)      
print(time.time()-start_time)    
