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
import h5py
import re
import sys

print(sys.argv)
exp=(sys.argv[1])
print(exp)
ri=int(sys.argv[2])
print(ri)

expdir = '/scratch/claire/speaker-listener/'
timeUnit='tr'
froidir='mor';

fnames=glob.glob(os.path.join(expdir + '/' + exp + '/fmri/timeseries/' + timeUnit + '/network/' + froidir , 'listenerAll_iscmasked_*.mat'))
rname=re.search(expdir + '/' + exp + '/fmri/timeseries/' + timeUnit + '/network/' + froidir + '/listenerAll_iscmasked_(.+?).mat', fnames[ri-1]).group(1) 

# load mat file
Ks_optimal = np.load(expdir +  exp+ '/fmri/hmm/Ks_optimal.npy',allow_pickle='TRUE').item()
K=Ks_optimal[rname]

data_mat=scipy.io.loadmat(expdir + '/' + exp + '/fmri/timeseries/' + timeUnit + '/network/' + froidir + '/zscore_speaker_'+ rname + '.mat')
data=data_mat['data']

data_mat=scipy.io.loadmat(expdir + '/' + exp + '/fmri/timeseries/' + timeUnit + '/network/' + froidir + '/zscore_listenerAll_' +rname+ '.mat')
gdata=data_mat['gdata']
g=np.mean(gdata,axis=2)

gdata_list=np.split(gdata.transpose((1, 0, 2)),gdata.shape[2],axis=2)
for i in range(len(gdata_list)):
    gdata_list[i]=np.squeeze(gdata_list[i])

tn=gdata.shape[1]
subjn=gdata.shape[2]
voxn=gdata.shape[0]

# Find the events in this dataset
seg = brainiak.eventseg.event.EventSegment(K)
seg.fit(g.T) 
# evi=seg.segments_[0].argmax(axis=1)+1
segments_LG=seg.segments_[0]
eventLabels_LG=seg.predict(g.T)

segmentsLG_S, _=seg.find_events(data.T);#, scramble=True)
eventLabelsLG_S=seg.predict(data.T)

segmentsLG_L= [None] * subjn
for si in range(subjn):
    segmentsLG_L[si],_=seg.find_events(gdata[:,:,si].T)

seg=seg.fit(gdata_list);#, scramble=True)
segments_L=seg.segments_

segmentsL_S, _=seg.find_events(data.T);
eventLabelsL_S=seg.predict(data.T)

segmentsL_S, _=seg.find_events(data.T);
eventLabelsL_S=seg.predict(data.T)

f = h5py.File(expdir +  exp+ '/fmri/hmm/'+ rname + '.hdf5', "w")
f.create_dataset('segments_LG',data=segments_LG)
f.create_dataset('segmentsLG_S',data=segmentsLG_S)
f.create_dataset('segmentsLG_L',data=segmentsLG_L)
f.create_dataset('segments_L',data=segments_L)
f.create_dataset('segmentsL_S',data=segmentsL_S)
f.create_dataset('eventLabels_LG',data=eventLabels_LG)
f.create_dataset('eventLabelsLG_S',data=eventLabelsLG_S)

f.close()
