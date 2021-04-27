function wholeBrain_temporal_lagcorr_SL_g_permPhase(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
crop_start=10;
lags_tested={-10:10, -30:30};
permN=1000;

% for ei=1:4;%1:4;
exp=experiments{ei};
mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/wholeBrain/SL_g/perm']);

for lagi=1;%1:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/speaker_zscore.mat' ],'data')
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll_zscore.mat' ],'gdata','keptvox');
    [voxn,tn,listenerN]=size(gdata);
    
            keptT=(crop_start+1):tn;
    
    y=nanmean(gdata(:,keptT,:),3);
    
    for perm=1:permN;
        rng(perm)
        data_perm=(phase_rand(data',0))';
        x=data_perm(:,keptT);
        
        r(:,:,perm)=(lagcorr_claire(y',x',lags))';
    end
    save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/wholeBrain/SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase'],'r','lags','keptvox','keptT','-v7.3');
end
% end




