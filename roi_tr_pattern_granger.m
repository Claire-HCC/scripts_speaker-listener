
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min
crop_start=10;
lags=-60:-4;
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
        Xs=reshape(Xs,roi_voxn*length(keptT),length(lags));
        
        [b_l(ri,:),~,r_l,~,stats_l]=regress(y,[ones(size(Xl,1),1) Xl]);
        [b_ls(ri,:),~,r_ls,~,stats_ls]=regress(y,[ones(size(Xl,1),1) Xl Xs]);
        [b_s(ri,:),~,r_s,~,stats_s]=regress(y,[ones(size(Xl,1),1) Xs]);
        
        r2_l(ri,:)=stats_l(1);
        F_l(ri,:)=stats_l(2);
        p_l(ri,:)=stats_l(3);
        
        r2_ls(ri,:)=stats_ls(1);
        F_ls(ri,:)=stats_ls(2);
        p_ls(ri,:)=stats_ls(3);
        
        
        r2_s(ri,:)=stats_s(1);
        F_s(ri,:)=stats_s(2);
        p_s(ri,:)=stats_s(3);
        
        %The numerator of the F-statistic
        rss_l = r_l'*r_l;
        rss_ls = r_ls'*r_ls;
        rss_s = r_s'*r_s;
        
        F_num = (rss_l- rss_ls)/length(lags);
        %The denominator of the F-statistic
        F_den = rss_ls/(length(y)-2*length(lags)-1);
        %The F-Statistic
        F_s2l(ri,1) = F_num/F_den;
        
        F_num = (rss_s- rss_ls)/length(lags);
        %The denominator of the F-statistic
        F_den = rss_ls/(length(y)-2*length(lags)-1);
        %The F-Statistic
        F_l2l(ri,1) = F_num/F_den;
        
        clear Xl Xs
    end
    
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/granger_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b_ls','b_l','b_s','F_s2l','F_l2l','F_ls','F_s','F_l','r2_l','r2_s','r2_ls','p_s','p_l','p_ls','lags','rnames','keptT');
    clear b F p r2 rnames coupling
end

toc
