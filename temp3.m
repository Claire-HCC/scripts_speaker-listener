
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

lags=-4:4;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');


for ei=2;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags))],'couplingz','rnames','keptT','b');
    cp_ll=nanmean(couplingz,3);
    keptT_ll=keptT;
    b_ll=b(:,2:end,:);
    
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags))],'b','couplingz','rnames','keptT','lags');
    b_sl=b(:,2:end);
    cp_sl=couplingz;
    keptT_sl=keptT;
    keptT=[max(min(keptT_ll),min(keptT_sl)):min(max(keptT_ll),max(keptT_sl))];
    
    herd_sig=zeros(length(rnames),1);
     load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm_nocrop' ],'couple_sig_fdr','couple_sig_fwe');
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm_nocrop' ],'couple_sig_fdr','couple_sig_fwe');
    
    for ri=1:length(rnames);
        [herd(ri,1) herd_p(ri,1)]=corr(cp_sl(ri,keptT)',cp_ll(ri,keptT)');
    end
    
    herd_sig(couple_sig_fdr==1,1)=fdr0(herd_p(couple_sig_fdr==1,1),0.05);
    
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(herd_sig(:,1)==1));
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, herd(herd_sig(:,1)==1,1)]);
    % save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr.nii']);
    
    % cluster betas by their temporal profiles instead of absolute value.
    b_sl_z=zscore(b_sl(:,:),0,2);
    % b_clusters(herd_sig(:,ei)==1,ei)=kmeans(b_sl_z(herd_sig(:,ei)==1,:),6,'Distance','correlation','Replicates',5);
    b_cluster(herd_sig(:,1)==1,1)=kmeans_opt(b_sl_z(herd_sig==1,:));
    
    
    nii=roiTable2wholeBrainNii_mor([roi_ids, b_cluster(herd_sig==1)]);
    % save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_b_clusters.nii']);
    
    rnames(herd_sig(:,1)==1)
    
    % save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_'  exp '_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'b_sl','b_ll','b_cluster','herd_sig','rnames','lags','herd','herd_p');
    
end
