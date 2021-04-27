import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, zscore, pearsonr, ttest_1samp
import os
import glob
import re

expdir = '/scratch/claire/speaker-listener/'
exp='sherlock'
timeUnit='tr'
froidir='mor'
fnames=glob.glob(os.path.join(expdir + '/' + exp + '/fmri/timeseries/' + timeUnit + '/network/' + froidir , 'listenerAll_iscmasked_*.mat'))
Ks_optimal={}

for fname in fnames:
    rname=re.search('/scratch/claire/speaker-listener//'+exp+'/fmri/timeseries/tr/network/mor/listenerAll_iscmasked_(.+?).mat', fname).group(1) 
    
    within_across_real=np.load(expdir +  exp+ '/fmri/hmm/WithinAcrossSim/'+ rname + '_ListenersLeave1Out_withinAcross.npy')   
    realzm=np.nanmean(np.arctanh(within_across_real),axis=1)
    Ks_optimal[rname] = np.nanargmax(realzm)+1
    # within_across_null=np.load(expdir +  exp+ '/fmri/hmm/WithinAcrossSim/'+ rname + '_ListenersLeave1Out_withinAcross_perm.npy')  
    # nullzm=np.nanmean(np.arctanh(within_across_null),axis=1)
    # t,_=ttest_1samp(nullzm,realzm,axis=1,nan_policy='omit')
    # real vs. null
    # t=-t
    # Ks_optimal[rname] = np.nanargmax(t)+1
    # plt.plot(t,'r')
    # plt.plot(np.nanmean(within_across_real,axis=1),'r')
    # plt.plot(np.nanmean(within_across_null,axis=1),'gray')
    # plt.title(rname +" "+str(np.nanargmax(np.nanmean(within_across_real,axis=1))+1))
    # plt.title(rname +" "+str(np.nanargmax(t)+1))
    # plt.show()
    # input("Press Enter to continue...")

np.save(expdir +  exp+ '/fmri/hmm/Ks_optimal.npy',Ks_optimal)

import csv
with open(expdir +  exp+ "/fmri/hmm/Ks_optimal.csv", 'w', newline="") as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in Ks_optimal.items():
        writer.writerow([key, value])
