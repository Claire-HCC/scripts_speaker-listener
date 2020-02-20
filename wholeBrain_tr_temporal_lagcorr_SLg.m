function wholeBrain_tr_temporal_lagcorr_SLg

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

crop_start=10;
lags_tested={-10:10, -30:30};

for ei=1:4;%1:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_speaker.mat' ],'data')
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_listenerAll.mat' ],'gdata','keptvox');
        [voxn,tn,listenerN]=size(gdata);
        
        keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
        keptT_e=min(tn,tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        x=data(:,keptT);
        y=nanmean(gdata(:,keptT,:),3);
        [r ]=lagcorr_claire(y',x',lags);
        r=r';
        
        save([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','keptvox','keptT','-v7.3');
    end
end




