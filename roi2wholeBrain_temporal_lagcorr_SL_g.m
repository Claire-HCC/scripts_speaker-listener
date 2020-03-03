function roi2wholeBrain_tr_temporal_lagcorr_SLg

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

crop_start=10;
lags_tested={-10:10, -30:30};
speakerSeed='vPCUN';

for ei=[1:2]
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' speakerSeed '.mat' ],'data')
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_listenerAll_grayMasked.mat' ],'gdata','keptvox');

        [voxn,tn,listenerN]=size(gdata);
        
        keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
        keptT_e=min(tn,tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        x=repmat(mean(data(:,keptT),1),length(keptvox),1);
        y=nanmean(gdata(:,keptT,:),3);
        
        [r]=lagcorr_claire(y',x',lags);
        r=r';
        clear y
    end

save([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi2wholeBrain/SLg/'  speakerSeed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','keptvox','keptT','-v7.3');
clear x y r
end


