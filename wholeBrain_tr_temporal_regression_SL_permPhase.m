function wholeBrain_tr_temporal_regression_SL_permPhase(perm)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

tic % 15 min

% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
lags=-10:10;

for ei=1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_listenerAll.mat' ],'gdata');
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_speaker.mat'],'data','keptvox');
    voxn=1;
    tn=size(gdata,2);
    
    data_perm=nan(size(data))';
    [data_perm((crop_start+1):end,:), ~]=phase_rand2(data(:,(crop_start+1):end)',1);
    data_perm=data_perm';
    
    keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
    keptT_e=min(tn,tn-max(lags));
    keptT=keptT_s:keptT_e;
    
    b=nan([voxn length(lags)+1 ]);
    r2=nan([voxn 1 ]);
    F=nan([voxn 1]);
    p=nan([voxn 1]);
    r=nan([voxn,tn]);
    
    for vi=1:size(data,1);
        g=nanmean(gdata(vi,:,:),3);
        y=g;
        y=y(:,keptT);
        y=y(:);
        
        for li=1:length(lags);
            X(:,:,li)=data_perm(vi,keptT+lags(li));
        end
        
        X=reshape(X,voxn*length(keptT),length(lags));
        
        % centralized X
        X=X-mean(X);
        
        % add an constant
        [b(vi,:),~,r(vi,keptT),~,stats]=regress(y,[ones(size(X,1),1) X]);
        
        r2(vi,1)=stats(1);
        F(vi,1)=stats(2);
        p(vi,1)=stats(3);

        clear X
    end
    
    % couplingz=0.5*log((1+coupling)./(1-coupling));
    save(sprintf('%s/%s/fmri/temporal_regression/%s/wholeBrain/SLg/perm/regression_SL_lag%d-%d_permPhase%04d',expdir,exp,timeUnit,min(lags),max(lags),perm),'b','F','r2','p','lags','keptT','r','-v7.3');
    
    clear b F p r2 r Y r2_byTime
end


