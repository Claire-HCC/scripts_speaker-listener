# coding: utf-8
# #!/bin/env python
import sys
import brainiak.eventseg.event
import numpy as np
import glob
import os
import scipy.io

# this should finish in 3hr for one roi

print(sys.argv)
exp=(sys.argv[1])
print(exp)
rname=(sys.argv[2])
print(rname)
K=int(sys.argv[3])
print(K)

expdir = '/scratch/claire/speaker-listener/'
timeUnit='tr'
froidir='mor';
exps=['pieman','bronx','merlin','sherlock']

fname=expdir+ exp + "/fmri/timeseries/tr/roi/mor/zscore_listenerAll_" + rname +".mat"
data_mat=scipy.io.loadmat(fname)
gdata=data_mat['gdata']

tn=gdata.shape[1]
subjn=gdata.shape[2]
voxn=gdata.shape[0]

eventLabels_others= np.empty((subjn,tn))
eventLabels_others[:]=np.nan
eventLabels_self= np.empty((subjn,tn))
eventLabels_self[:]=np.nan

for s in range(subjn):
    if np.sum(~np.isnan(gdata[:,:,s]))>0:
        othersi=np.arange(subjn)
        othersi=othersi[othersi!=s]
        others=np.nanmean(gdata[:,:,othersi],axis=2)
        self=gdata[:,:,s]
        
        # Find the events in this dataset
        seg = brainiak.eventseg.event.EventSegment(K)
        seg.fit(others.T) 
        eventLabels_others[s,:]=seg.predict(others.T)+1
        
        segments, _=seg.find_events(self.T);#, scramble=True)
        eventLabels_self[s,:]=np.nanargmax(segments,axis=1)+1
    print(s)

np.savez(expdir +  exp+ '/fmri/hmm/'+ rname + '_ListenersLeave1Out_eventLabels_K' +str(K)+ '.npz', eventLabels_others=eventLabels_others,eventLabels_self=eventLabels_self) 


