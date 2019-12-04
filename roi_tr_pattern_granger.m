
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min
crop_start=10;
lags=-60:-5;
for ei=3;%1:2;
    exp=experiments{ei};
    rnames=dir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/'  froidir '/zscore_listenerAll_*.mat']);
    rnames={rnames.name};
    rnames=strrep(rnames,'zscore_listenerAll_','');
    rnames=strrep(rnames,'.mat','');
    rnames=rnames';
    
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
        
        g=nanmean(gdata,3);
        y=g;
        
        keptT_s=find(([1:roi_tn]+min(lags))==1)++crop_start;
        keptT_e=min(roi_tn,roi_tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        y=y(:,keptT);
        y=y(:);
        
        for li=1:length(lags);
            Xs(:,:,li)=data(:,keptT+lags(li));
            Xl(:,:,li)=g(:,keptT+lags(li));
        end
        
        Xl=reshape(Xl,roi_voxn*length(keptT),length(lags));
        Xl=[ones(size(Xl,1),1) Xl];
        
        Xs=reshape(Xs,roi_voxn*length(keptT),length(lags));
       
        [b_l(ri,:),~,~,~,stats_l]=regress(y,Xl);
        [b_ls(ri,:),~,~,~,stats_ls]=regress(y,[Xl Xs]);
        
        r2_l(ri,:)=stats_l(1);
        F_l(ri,:)=stats_l(2);
        p_l(ri,:)=stats_l(3);
        
        r2_ls(ri,:)=stats_ls(1);
        F_ls(ri,:)=stats_ls(2);
        p_ls(ri,:)=stats_ls(3);
        
        clear Xl Xs
    end
    
      save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/granger_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b_l','b_ls','F_l','r2_l','p_ls','F_ls','r2_ls','p_ls','lags','rnames','keptT');
    clear b F p r2 rnames coupling
end

toc
