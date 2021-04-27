function wholeBrain_temporal_regression_SL_g(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

crop_start=10;
lags_tested={-4:4, -10:10};


exp=experiments{ei};
% mkdir([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_g/perm'])
for lagi=1:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/speaker_zscore.mat' ],'data')
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll_zscore.mat' ],'gdata','keptvox');
    gdata(:,:,subjects_excluded{ei})=NaN;
    [voxn,tn,listenerN]=size(gdata);
    
    keptT_s= max(-min(lags),0)+crop_start+1;
    keptT_e=min(tn,tn-max(lags));
    keptT=keptT_s:keptT_e;
    
    b=nan([length(keptvox) length(lags) ]);
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
        
        [btemp,~,r(vi,keptT),~,stats]=regress(y,[ones(size(X,1),1) X]);
        b(vi,:)=btemp(2:end);
        r2(vi,1)=stats(1);
        F(vi,1)=stats(2);
        p(vi,1)=stats(3);
        
        clear X
    end
    
    save([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','keptT','keptvox','-v7.3');
    clear b F p r2 r2_byTime
end


