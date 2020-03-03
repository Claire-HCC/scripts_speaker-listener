function roi2wholeBrain_tr_temporal_lagcorr_SLg_permPhase(perm)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

crop_start=10;
lags_tested={-10:10, -30:30};
speakerSeed='vPCUN';


for ei=[1:4];
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load(sprintf('%s/%s/fmri/timeseries/%s/roi/%s/permPhase/zscore_speaker_%s.mat',expdir,exp,timeUnit,froidir,speakerSeed),'data');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_listenerAll_grayMasked.mat' ],'gdata','keptvox');

        [voxn,tn,listenerN]=size(gdata);
        
        keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
        keptT_e=min(tn,tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        x=repmat(mean(data(:,keptT,perm),1),length(keptvox),1);
        y=nanmean(gdata(:,keptT,:),3);
        
        [r]=lagcorr_claire(y',x',lags);
        r=r';
        clear y
    end
    
    save(sprintf('%s/%s/fmri/temporal_lagcorr/%s/roi2wholeBrain/SLg/perm/%s_lag%d-%d_permPhase%04d',expdir,exp,timeUnit,speakerSeed,min(lags),max(lags),perm),'r','lags','keptvox','keptT','keptvox','-v7.3');
    clear x y r
end


