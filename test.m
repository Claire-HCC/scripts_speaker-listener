
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';

exp='merlin';
load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_dPCC' ],'data');
crop_start=10;

roi_voxn=size(data,1);
roi_tn=size(data,2);

lags=-20:-1;

b=[];
r=[];

for vi=1:size(data,1);
    keptT=(abs(min(lags))+1):(roi_voxn-max(lags));
    % substract the global mean pattern
    y=data(vi,keptT)';
    
    for li=1:length(lags);
        X(:,li)=data(vi,keptT+lags(li));
    end
    
    % centralized X
    X=X-mean(X);
    
    % add an constant
    X=[ones(size(X,1),1) X];
    
    [b(vi,:),~,r(vi,:)]=regress(y,X);
    clear X y
end
