function wholeBrain_temporal_regression_SL_g_permPhase(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
lags_tested={-4:4, -10:10};
permN=1000;%

exp=experiments{ei};

mkdir(sprintf('%s/%s/fmri/temporal/regression/%s/wholeBrain/SL_g/perm/',expdir,exp,timeUnit))
for lagi=2;%1:length(lags_tested);
    lags=lags_tested{lagi};
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll_zscore.mat' ],'gdata');
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/speaker_zscore.mat'],'data','keptvox');
    gdata(:,:,subjects_excluded{ei})=NaN;
    g=nanmean(gdata,3);
    voxn=1;
    tn=size(gdata,2);
    
    keptT_s= max(-min(lags),0)+crop_start+1;
    keptT_e=min(tn,tn-max(lags));
    keptT=keptT_s:keptT_e;
    
    b=nan([voxn length(lags)+1 permN]);
    r2=nan([voxn permN ]);
    F=nan([voxn permN]);
    p=nan([voxn permN]);
    
    for perm=1:permN;
        rng(perm);
        
        data_perm=nan(size(data))';
        [data_perm((crop_start+1):end,:)]=phase_rand2(data(:,(crop_start+1):end)',1);
        data_perm=data_perm';
        
        for vi=1:size(data,1);
            
            y=g(vi,keptT);
            y=y(:);
            
            for li=1:length(lags);
                X(:,:,li)=data_perm(vi,keptT+lags(li));
            end
            
            X=reshape(X,voxn*length(keptT),length(lags));
            
            % centralized X
            X=X-mean(X);
            
            % add an constant
            [b(vi,:,perm),~,~,~,stats]=regress(y,[ones(size(X,1),1) X]);
            
            r2(vi,perm)=stats(1);
            F(vi,perm)=stats(2);
            p(vi,perm)=stats(3);
            
            clear X
        end
    end

    save(sprintf('%s/%s/fmri/temporal/regression/%s/wholeBrain/SL_g/perm/lag%d-%d_permPhase',expdir,exp,timeUnit,min(lags),max(lags)),'b','F','r2','p','lags','keptT','-v7.3');
    
    clear b F p r2 Y r2_byTime
end




