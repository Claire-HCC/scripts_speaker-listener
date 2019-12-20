clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
crop_start=10;
lags=[-10:-4];

for ei=[3];%1:2;
    exp=experiments{ei};
    
    rnames=table2array(roi_table(:,3));
    ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
    rnames=rnames(ris);
    
    for ri=1:length(rnames);
        clear data_mat
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        
        keptT_s=find(([1:roi_tn]+min(lags))==1)+crop_start;
        keptT_e=min(roi_tn,roi_tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        g=nanmean(gdata,3);
        y=g;
        y=y(:,keptT);
        y=y(:);
        
        for li=1:length(lags);
            % lagmatrix
            Xs(:,:,li)=data(:,keptT+lags(li));
            Xl(:,:,li)=g(:,keptT+lags(li));
        end
        
        
        Xl=reshape(Xl,roi_voxn*length(keptT),length(lags));
        Xs=reshape(Xs,roi_voxn*length(keptT),length(lags));
        
        % centralized X
        Xl=Xl-mean(Xl);
        Xs=Xs-mean(Xs);
        
        
        [b_l(ri,1:(1+length(lags))),~,r_l,~,stats_l]=regress(y,[ones(size(Xl,1),1) Xl]);
        [b_ls(ri,1:(1+length(lags)*2)),~,r_ls,~,stats_ls]=regress(y,[ones(size(Xl,1),1) Xl Xs]);
        [b_s(ri,1:(1+length(lags))),~,r_s,~,stats_s]=regress(y,[ones(size(Xl,1),1) Xs]);
        
        r2_l(ri,1)=stats_l(1);
        F_l(ri,1)=stats_l(2);
        p_l(ri,1)=stats_l(3);
        
        r2_ls(ri,1)=stats_ls(1);
        F_ls(ri,1)=stats_ls(2);
        p_ls(ri,1)=stats_ls(3);
        
        r2_s(ri,1)=stats_s(1);
        F_s(ri,1)=stats_s(2);
        p_s(ri,1)=stats_s(3);
        
        
        %The numerator of the F-statistic
        rss_l = r_l'*r_l;
        rss_s = r_s'*r_s;
        rss_ls = r_ls'*r_ls;
        
        F_num = (rss_l- rss_ls)/length(lags);
        %The denominator of the F-statistic
        F_den = rss_ls/(length(y)-2*length(lags)-1);
        %The F-Statistic
        F_s2l(ri,1) = F_num/F_den;
        
        F_num = (rss_s- rss_ls)/length(lags);
        F_l2l(ri,1) = F_num/F_den;
        
        
        ymat=reshape(y,roi_voxn,length(keptT));
        rmat_l=reshape(r_l,roi_voxn,length(keptT));
        rmat_ls=reshape(r_ls,roi_voxn,length(keptT));
        rmat_s=reshape(r_s,roi_voxn,length(keptT));
        ssr_l=sum(rmat_l.^2);
        ssr_ls=sum(rmat_ls.^2);
        ssr_s=sum(rmat_s.^2);
        sst=sum((ymat-mean(y)).^2);
        r2_l_byTime(ri,:)=nan(roi_tn,1);
        r2_l_byTime(ri,keptT)=1-ssr_l./sst;
        r2_ls_byTime(ri,:)=nan(roi_tn,1);
        r2_ls_byTime(ri,keptT)=1-ssr_ls./sst;
        r2_s_byTime(ri,:)=nan(roi_tn,1);
        r2_s_byTime(ri,keptT)=1-ssr_s./sst;
        
        clear Xl Xs
    end
    
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/granger_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ''],'F_s2l','lags','rnames','keptT',...
        'r2_l_byTime','r2_ls_byTime','r2_s_byTime');
    clear F_s2l
end



toc
