function wholeBrain_tr_temporal_lagcorr_SLg_permPhase(perm)

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
        
        load( sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/perm/zscore_speaker_permPhase%04d.mat',expdir, exp,timeUnit,perm),'data')
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_listenerAll.mat' ],'gdata','keptvox');
        [voxn,tn,listenerN]=size(gdata);
        
        keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
        keptT_e=min(tn,tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        x=data(:,keptT);
        y=nanmean(gdata(:,keptT,:),3);
        [r ]=lagcorr_claire(y',x',lags);
        r=r';
        save(sprintf('%s/%s/fmri/temporal_lagcorr/%s/wholeBrain/SLg/perm/lag%d-%d_permPhase%04d',expdir,exp,timeUnit,min(lags),max(lags),perm),'r','lags','keptvox','keptT','-v7.3');
    end
    
end



