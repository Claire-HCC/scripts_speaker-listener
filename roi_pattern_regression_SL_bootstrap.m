clear all;

tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
lags=-4:4;
bootsN=1000;

for ei=3;
    exp=experiments{ei};
    rnames=dir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/'  froidir '/zscore_listenerAll_*.mat']);
    rnames={rnames.name};
    rnames=strrep(rnames,'zscore_listenerAll_','');
    rnames=strrep(rnames,'.mat','');
    rnames=rnames';
    rnames=rnames(ismember(rnames,roi_table.region));
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rnames{1} '.mat' ],'gdata');
    g1=combnk(1:size(gdata,3),9);
    g1=g1(randperm(size(g1,1)),:);
    
    b=[];
    roi_ids=[];
    r2=[];
    F=[];
    for  booti=1:bootsN;
        for ri=1:length(rnames);
            clear data_mat
            rname=rnames{ri};
            
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
            
            roi_voxn=size(gdata,1);
            roi_tn=size(gdata,2);
            
            g=nanmean(gdata(:,:,g1(booti,:)),3);
            y=g;
            keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags, 0]));
            y=y(:,keptT);
            y=y(:);
            
            for li=1:length(lags);
                X(:,:,li)=data(:,keptT+lags(li));
            end
            
            X=reshape(X,roi_voxn*length(keptT),length(lags));
            
            % centralized X
            X=X-mean(X);
            
            % add an constant
            X=[ones(size(X,1),1) X];
            
            [b(ri,:,1,booti),~,~,~,stats]=regress(y,X);
            
            r2(ri,:,1,booti)=stats(1);
            F(ri,:,1,booti)=stats(2);
            p(ri,:,1,booti)=stats(3);
            
            Y=X*b(ri,:,1,booti)';
            Y=reshape(Y,roi_voxn,length(keptT));
            y=reshape(y,roi_voxn,length(keptT));
            coupling(ri,:,1,booti)=zeros(1,roi_tn);
            coupling(ri,keptT,:,booti)=corr_col(y,Y);
            
            clear X
        end
    end
    couplingz=0.5*log((1+coupling)./(1-coupling));
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_bootstrap' ],'b','F','r2','p','couplingz','lags','rnames','keptT','g1','bootsN');
    clear b F p r2 rnames coupling
end

toc
beep
