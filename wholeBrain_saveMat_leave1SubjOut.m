close all
clear all
loc='cluster';

set_parameters

tic

for ei=2;%1:2;
    exp=experiments{ei};
    load([expdir '/' exp '/' '/fmri/timeseries/tr/wholeBrain/zscore_listenerAll.mat'],'gdata','keptvox');
    subjN=size(gdata,3);
    gdata_self=gdata;
    clear gdata
    
    for si=1:subjN;
        othersi=1:subjN;
        othersi=othersi(~ismember(othersi,si));
        
        gdata(:,:,si)=mean(gdata_self(:,:,othersi),3);
        
    end
    save([expdir '/' exp '/' '/fmri/timeseries/tr/wholeBrain/zscore_listenerAll_others.mat'],'gdata','keptvox','-v7.3');
end