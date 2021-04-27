function wholeBrain_temporal_regression_SL_each(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

crop_start=10;
lags_tested={-4:4, -10:10};

exp=experiments{ei};
% mkdir([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_each/perm'])
for lagi=2;%1:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/speaker_zscore.mat' ],'data')
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll_zscore.mat' ],'gdata','keptvox');
    % gdata(:,:,subjects_excluded{ei})=NaN;
    [voxn,tn,listenerN]=size(gdata);
    
    keptT_s= max(-min(lags),0)+crop_start+1;
    keptT_e=min(tn,tn-max(lags));
    keptT=keptT_s:keptT_e;
    
    b=nan([length(keptvox) length(lags) listenerN ]);
    r2=nan([length(keptvox) 1 listenerN]);
    F=nan([length(keptvox) 1 listenerN]);
    p=nan([length(keptvox) 1 listenerN]);
    r2_byTime=nan([length(keptvox),tn, listenerN]);
    
    for vi=1:length(keptvox);
        
        for li=1:length(lags);
            X(:,:,li)=data(vi,keptT+lags(li));
        end
        X=reshape(X,length(keptT),length(lags));
        X=X-mean(X);
        
        for si=1:listenerN;
            y=gdata(vi,keptT,si)';
            
            [btemp,~,r(vi,keptT,si),~,stats]=regress(y,[ones(size(X,1),1) X]);
            b(vi,:,si)=btemp(2:end);
            r2(vi,1,si)=stats(1);
            F(vi,1,si)=stats(2);
            p(vi,1,si)=stats(3);
            
        end
        clear X
    end
    
    save([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','keptT','keptvox','-v7.3');
    clear b F p r2 r2_byTime
end


