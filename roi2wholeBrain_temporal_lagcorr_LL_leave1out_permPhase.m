function roi2wholeBrain_tr_temporal_lagcorr_LL_leave1out_permPhase(ei);

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
crop_start=10;
lags_tested={-10:10, -30:30};
seed='vPCUN';
permN=10000;

exp=experiments{ei};

for lagi=1;%1:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/perm/listenerAll_zscore_' seed '.mat' ],'gdata','keptvox');
    gdata_seed=mean(gdata,1);
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll_zscore.mat' ],'gdata','keptvox');
    [voxn,tn,listenerN]=size(gdata);
    
    keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
    keptT_e=min(tn,tn-max(lags));
    keptT=keptT_s:keptT_e;
    
    r=nan([voxn length(lags) listenerN permN]);
    for si=1:listenerN;
        othersi=1:listenerN;
        othersi=othersi(othersi~=si);
        y=nanmean(gdata(:,keptT,othersi),3);
        
        for perm=permN;
            x=repmat(gdata_seed(:,keptT,si,perm),length(keptvox),1);
            r(:,:,si,perm)=(lagcorr_claire(y',x',lags))';
        end
    end
    
    save(sprintf('%s/%s/fmri/temporal/lagcorr/%s/roi2wholeBrain/LL_leave1out/perm/%s_lag%d-%d_permPhase',expdir,exp,timeUnit,speakerSeed,min(lags),max(lags)),'r','lags','keptT','keptvox','-v7.3');
    clear r
end



