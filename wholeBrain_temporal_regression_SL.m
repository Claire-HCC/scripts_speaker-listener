function wholeBrain_tr_temporal_regression_SL(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

crop_start=10;
lags_tested={-10:10};


exp=experiments{ei};

for lagi=1%:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_speaker.mat' ],'data')
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_listenerAll.mat' ],'gdata','keptvox');
    [voxn,tn,listenerN]=size(gdata);
    
    keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
    keptT_e=min(tn,tn-max(lags));
    keptT=keptT_s:keptT_e;
    
    b=nan([length(keptvox) length(lags)+1 ]);
    r2=nan([length(keptvox) 1 ]);
    F=nan([length(keptvox) 1]);
    p=nan([length(keptvox) 1]);
    r2_byTime=nan([length(keptvox),tn]);
    
    g=nanmean(gdata,3);
    
    for vi=1:length(keptvox);
        y=g(vi,keptT)';
        
        for li=1:length(lags);
            X(:,:,li)=data(vi,keptT+lags(li));
        end
        
        X=reshape(X,length(keptT),length(lags));
        X=X-mean(X);
        
        [b(vi,:),~,r(vi,keptT),~,stats]=regress(y,[ones(size(X,1),1) X]);
        
        r2(vi,1)=stats(1);
        F(vi,1)=stats(2);
        p(vi,1)=stats(3);
        
        clear X
    end
    
    save([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','keptT','keptvox','-v7.3');
    clear b F p r2 r2_byTime
end


