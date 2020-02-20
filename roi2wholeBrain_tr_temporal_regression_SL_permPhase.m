function roi2wholeBrain_tr_temporal_regression_SL_permPhase(perm)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
speakerSeed='vPCUN';
crop_start=10;
lags_tested={-10:10};

for ei=1:4;%2:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load(sprintf('%s/%s/fmri/timeseries/%s/roi/%s/permPhase/speaker_%s.mat',expdir,exp,timeUnit,froidir,speakerSeed),'data');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll.mat' ],'gdata','keptvox');
        [voxn,tn,listenerN]=size(gdata);
        
        b=nan([length(keptvox) length(lags)+1 ]);
        r2=nan([length(keptvox) 1 ]);
        F=nan([length(keptvox) 1]);
        p=nan([length(keptvox) 1]);
        r=nan(size(data));
        
        keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
        keptT_e=min(tn,tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        g=nanmean(gdata,3);
        data=nanmean(data(:,:,perm),1);
        
        for vi=1:length(keptvox);
            y=g(vi,keptT)';
            
            for li=1:length(lags);
                X(:,:,li)=data(1,keptT+lags(li));
            end
            
            X=reshape(X,length(keptT),length(lags));
            X=X-mean(X);
            
            [b(vi,:),~,r(vi,keptT),~,stats]=regress(y,[ones(size(X,1),1) X]);
            
            r2(vi,1)=stats(1);
            F(vi,1)=stats(2);
            p(vi,1)=stats(3);
            
            clear X
        end
       
        save(sprintf('%s/%s/fmri/temporal_regression/%s/roi2wholeBrain/%s/SLg/perm/regression_SL_lag%d-%d_permPhase%04d',expdir,exp, timeUnit,froidir,min(lags),max(lags),perm),'b','F','r2','p','lags','keptT','-v7.3');
        clear b F p r2 r
    end
end

