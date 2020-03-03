clear all;

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
crop_start=10;
iters=1000;
lags=0;

for ei=3;
    exp=experiments{ei};
    exp=experiments{ei};
    rnames=table2array(roi_table(:,3));
    ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
    rnames=rnames(ris);
    
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','rnames','keptT','r2_byTime');
    b_real=mean(b(:,2:end,:),3);
    r2_real=mean(r2,2);
    r2_byTime_real=mean(r2_byTime,3);
    tn=size(r2_byTime,2);
    
    for iter=1:iters;
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' num2str(iter)],'b','F','r2','p','lags','rnames','keptT','r2_byTime');
        b_null(:,:,iter)=mean(b(:,2:end,:),3);
        r2_byTime_null(:,:,iter)=mean(r2_byTime,3);
        r2_null(:,iter)=mean(r2,2);
    end
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        [~,r2_p(ri,1),~,stats]=ttest(r2(ri,:)',0,'tail','right');
        r2_t(ri,1)=stats.tstat;
    end
    r2_sig_fdr=(fdr0(r2_p,0.05)==1);
    r2_sig_fwe=(r2_p<(0.05/length(p))==1);
    
    r2_byTime_p=nan(length(rnames),tn);
    t_r2_byTime= r2_byTime_p;
    for ri=1:length(rnames);
        r2_byTime_p=nan(size(r2_byTime,1),size(r2_byTime,2));
        [~,r2_byTime_p(ri,:,:),~,stats] = ttest(squeeze(r2_byTime(ri,:,:))',0,'tail','right');
        r2_byTime_t(ri,:,:)=stats.tstat;
    end
    i=zeros(size(r2_byTime_p));
    i(r2_sig_fdr==1,keptT)=1;
    i=i(:);
    r2_byTime_sig_fdr=nan(size(r2_byTime_p));
    
    for ri=1:length(rnames);
        if r2_sig_fdr(ri)==1;
            r2_byTime_sig_fdr(ri,keptT)=(fdr0(r2_byTime_p(ri,keptT),0.05)==1);
        end
    end
    
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_statsPhase.mat'],'b_real','b_null','b_p','b_t','b_sig_fdr','r2_p',...
        'r2_t','r2_sig_fdr','r2_byTime_p','r2_byTime_t','r2_byTime_sig_fdr','lags','rnames','keptT','r2_byTime_null','r2_byTime_real');
    
    clear b_null b_real beta_p F_real F_null b_t b_p
end


