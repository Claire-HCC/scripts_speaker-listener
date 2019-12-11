clear all;
% F-test: https://support.sas.com/rnd/app/ets/examples/granger/index.htm
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

for ei=3;%1:4;
    exp=experiments{ei};
    
    rnames=dir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/'  froidir '/zscore_listenerAll_*.mat']);
    rnames={rnames.name};
    rnames=strrep(rnames,'zscore_listenerAll_','');
    rnames=strrep(rnames,'.mat','');
    rnames=rnames';
    rnames=rnames(ismember(rnames,roi_table.region));
    
    for ri=1:length(rnames);
      
        rname=rnames{ri};
        eventLabels_LG(ri,:)=h5read([expdir   exp '/fmri/hmm/findListenersEventInSpeaker/' rname  '.hdf5'],'///eventLabels_LG');
      
    end
end