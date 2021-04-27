clear all;

loc='cluster';
set_parameters;

for ei =1;
    exp=experiments{ei};
    tic % 15 min
    mkdir([expdir '/'  exp '/fmri/timeseries/tr/wholeBrain/']);
    load([expdir '/'  exp '/fmri/timeseries/tr_uncropped/wholeBrain/listenerAll.mat'],'gdata','fs','keptvox');

    crop_start=crop_start_listener{ei};
    voln=voln_listener{ei};
    gdata=gdata(:,(crop_start+1):(crop_start+voln),:);
    
    save([expdir '/'  exp '/fmri/timeseries/tr/wholeBrain/listenerAll.mat'],'gdata','fs','keptvox','-v7.3');
    
   % gdata=zscore(gdata,0,2);
   % save([expdir '/'  exp '/fmri/timeseries/tr/wholeBrain/listenerAll_zscore.mat'],'gdata','fs','keptvox','-v7.3');
    
    clear gdata
end
toc
