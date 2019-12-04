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
lags=-60:60;
for ei=3;
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
    for ri=1:length(rnames);
        clear data_mat
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');

        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        
        keptT_s=find(([1:roi_tn]+min(lags))==1)++crop_start;
        keptT_e=min(roi_tn,roi_tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        g=nanmean(gdata,3);
        y=g;
        y=y(:,keptT);
        y=y(:);
        
        for li=1:length(lags);
            X(:,:,li)=data(:,keptT+lags(li));
        end
        
        X=reshape(X,roi_voxn*length(keptT),length(lags));
        X=[ones(size(X,1),1) X];

        [b(ri,:),~,~,~,stats]=regress(y,X);
        
        r2(ri,:)=stats(1);
        F(ri,:)=stats(2);
        p(ri,:)=stats(3);
        
        Y=X*b(ri,:)';
        Y=reshape(Y,roi_voxn,length(keptT));
        
        clear X
    end
    
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','rnames','keptT');
    clear b F p r2 rnames coupling
end

toc
