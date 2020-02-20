function roi2wholeBrain_tr_temporal_lagcorr_SLg_permPhase(perm)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

crop_start=10;
lags_tested={-10:10, -30:30};
speakerSeed='vPCUN';

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load(sprintf('%s/%s/fmri/timeseries/%s/roi/%s/permPhase/speaker_%s.mat',expdir,exp,timeUnit,froidir,speakerSeed),'data');
        data=data(:,:,perm);
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll.mat' ],'gdata','keptvox');
        [voxn,tn,listenerN]=size(gdata);
        g=nanmean(gdata,3);
        
        keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
        keptT_e=min(tn,tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        x=mean(data(:,keptT),1);
        r=nan([voxn  length(lags) ]);
        
        for vi=1:voxn;
            
            y=g(vi,keptT);
            [r(vi,:) ,~]=lagcorr(y',x',lags);
            
            clear y
        end
    end
    mkdir([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi2wholeBrain/' froidir '/SLg/perm']);
    
    save(sprintf('%s/%s/fmri/temporal_lagcorr/%s/roi2wholeBrain/%s/SLg/perm/%s_lag%d-%d_permPhase%04d',expdir,exp,timeUnit,froidir,speakerSeed,min(lags),max(lags),perm),'r','lags','keptvox','keptT','keptvox','-v7.3');
    clear x y r
end


