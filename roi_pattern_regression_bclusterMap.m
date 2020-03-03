clear all;
close all
loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
lags=-4:4;
smoothSpan=3;

for ei=1:4;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F','b','rnames');
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],'F_perm','b_perm','rnames','couple_sig_fdr');
    
    bz=zscore(b(:,2:end),0,2);
    b_clusters=zeros(length(rnames),1);
    b_clusters(couple_sig_fdr==1)=kmeans_opt(bz(couple_sig_fdr==1,2:end));
        
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(couple_sig_fdr==1));
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, b_clusters(couple_sig_fdr==1)]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_b_clusters.nii']);
    
    
    for ci=1:max(b_clusters);
        roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(b_clusters==ci));
        roi_ids=cell2mat(roi_table.id(roi_table_inds));
        nii=roiTable2wholeBrainNii_mor([roi_ids, ci*ones(length(roi_ids),1)]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_cluster' num2str(ci) '.nii']);
        save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],'b_clusters','couple_sig_fdr','couple_sig_fwe','F_perm','b_perm','rnames');
    end
   
    for ci=1:max(b_clusters);
        figure;
        plot(lags*tr(ei),mean(bz(couple_sig_fdr==1 & b_clusters==ci,:),1),'k','linewidth',2);
        xlabel(['speaker precedes      shift (sec)    listener precedes']);
        ylabel('mean weigths');
        set(gca,'fontsize',14)
        grid on
        print(gcf,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_cluster' num2str(ci) '.tif'],'-dtiff');
        
    end
    
end
