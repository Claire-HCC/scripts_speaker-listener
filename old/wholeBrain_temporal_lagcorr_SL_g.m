function wholeBrain_temporal_lagcorr_SL_g

loc='cluster';
set_parameters;
timeUnit='tr' ;
crop_start=10;
lags_tested={-4:4, -10:10};

for ei=1:4;%1:4;
    exp=experiments{ei};
    mkdir([[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/wholeBrain/SL_g/perm']);
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/speaker_zscore.mat' ],'data')
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll_zscore.mat' ],'gdata','keptvox');
        gdata(:,:,subjects_excluded{ei})=NaN;
        [voxn,tn,listenerN]=size(gdata);
        
        keptT=(crop_start+1):tn;
        
        y=nanmean(gdata(:,keptT,:),3);
        
        x=data(:,keptT);
        r=(lagcorr_claire(y',x',lags))';
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','keptvox','keptT','-v7.3');
    end
end




