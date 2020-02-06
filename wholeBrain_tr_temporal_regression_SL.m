function wholeBrain_tr_temporal_regression_SL(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

tic % 15 min

% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
lags=-10:10;


exp=experiments{ei};

mkdir([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/wholeBrain/' froidir '/SLg/perm/']);

load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_listenerAll.mat' ],'gdata');
load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_speaker.mat'],'data','keptvox');
voxn=1;
tn=size(gdata,2);

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
        X(:,:,li)=data(vi,keptT+lags(li));
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
save([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/wholeBrain/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','keptT','r','keptvox','-v7.3');
clear b F p r2 r Y r2_byTime


