# coding: utf-8
# #!/bin/env python
import sys
import brainiak.eventseg.event
import numpy as np
import glob
import os
import scipy.io
from scipy.stats import stats
import logging
import matplotlib.pyplot as plt
from scipy.stats import norm, zscore, pearsonr
import pickle
import time
# this should finish in 3hr for one roi

rname=(sys.argv[1])
exp=(sys.argv[2])
expdir = '/scratch/claire/speaker-listener/'
timeUnit='tr'
froidir='mor';
exps=['pieman','bronx','merlin','sherlock']
eventN_test=range(10,121)

for fname in  glob.glob(os.path.join(expdir + '/' + exp + '/fmri/timeseries/' + timeUnit + '/roi/' + froidir , 'zscore_listenerAll_'+rname+'.mat')):
    data_mat=scipy.io.loadmat(fname)
    gdata=data_mat['gdata']
    
    tn=gdata.shape[1]
    subjn=gdata.shape[2]
    voxn=gdata.shape[0]
    
    eventLabels_others= np.empty((120,tn,subjn))
    eventLabels_others[:]=np.nan
    eventLabels_self= np.empty((120,tn,subjn))
    eventLabels_self[:]=np.nan
    
    for s in range(subjn):
        if np.sum(~np.isnan(gdata[:,:,s]))>0:
            othersi=np.arange(subjn)
            othersi=othersi[othersi!=s]
            others=np.nanmean(gdata[:,:,othersi],axis=2)
            self=gdata[:,:,s]
             
            for K in eventN_test:
                Ki=K-1
                
                # Find the events in this dataset
                seg = brainiak.eventseg.event.EventSegment(K)
                seg.fit(others.T) 
                eventLabels_others[Ki,:,s]=seg.predict(others.T)+1
                
                segments, _=seg.find_events(self.T);#, scramble=True)
                eventLabels_self[Ki,:,s]=np.nanargmax(segments,axis=1)+1
        print(s)
        
    np.savez(expdir +  exp+ '/fmri/hmm/'+ rname + '_ListenersLeave1Out_eventLabels.npz', eventLabels_others=eventLabels_others,eventLabels_self=eventLabels_self) 


