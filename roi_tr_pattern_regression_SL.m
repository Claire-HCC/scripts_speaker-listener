clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;

lags=-10:10;
for ei=1:4;
    exp=experiments{ei};
    rnames=dir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/'  froidir '/zscore_listenerAll_*.mat']);
    rnames={rnames.name};
    rnames=strrep(rnames,'zscore_listenerAll_','');
    rnames=strrep(rnames,'.mat','');
    rnames=rnames';
    rnames=rnames(ismember(rnames,roi_table.region));
    
    b=[];
    roi_ids=[];
    r2=[];
    F=[];
    for ri=54;%1:length(rnames);
        clear data_mat
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
        
        %% isc or pattern
        % gdata=nanmean(gdata,1);
        % data=nanmean(data,1);
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        
        g=nanmean(gdata,3);
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
        
        [b(ri,:),~,~,~,stats]=regress(y,X);
        
        r2(ri,:)=stats(1);
        F(ri,:)=stats(2);
        p(ri,:)=stats(3);
        
        Y=X*b(ri,:)';
        Y=reshape(Y,roi_voxn,length(keptT));
        y=reshape(y,roi_voxn,length(keptT));
        coupling(ri,:)=zeros(1,roi_tn);
        coupling(ri,keptT)=corr_col(y,Y);
        
        clear X
    end
    
    couplingz=0.5*log((1+coupling)./(1-coupling));
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','couplingz','lags','rnames','keptT');
    clear b F p r2 rnames coupling
end

toc
beep
