clear all;
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

binSize=30;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
lags=-4:4;

for ei=3%:4;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/bined_pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_bin' num2str(binSize) '_lag'  num2str(min(lags)) '-' num2str(max(lags)) ],'F','b','rnames');
    load([expdir '/' exp '/fmri/bined_pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_bin' num2str(binSize) '_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],'F_perm','b_perm','rnames');
    
    
    p=sum(F<F_perm,3)/1000;
    couple_sig_fwe(:,1)=(p<(0.05/length(rnames)));
    couple_sig_fdr(:,1)=fdr0(p,0.05);
    couple_sig_fwemax=(F>quantile(max(squeeze(F_perm)),0.95));

    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(couple_sig_fwe(:,1)==1));
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, F(couple_sig_fwe(:,1)==1,1)]);
    save_nii(nii,[expdir '/' exp '/fmri/bined_pattern_regression/' timeUnit '/roi/' froidir '/coulpingF_sl_bin' num2str(binSize) 'lag' num2str(min(lags)) '-' num2str(max(lags)) '_fwe.nii']);
    
    rnames(couple_sig_fdr(:,1)==1)
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(couple_sig_fdr(:,1)==1));
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, F(couple_sig_fdr(:,1)==1,1)]);
    save_nii(nii,[expdir '/' exp '/fmri/bined_pattern_regression/' timeUnit '/roi/' froidir '/coulpingF_sl_bin' num2str(binSize) 'lag' num2str(min(lags)) '-' num2str(max(lags)) '_fdr.nii']);
    
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(couple_sig_fwemax(:,1)==1));
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, F(couple_sig_fwemax(:,1)==1,1)]);
    save_nii(nii,[expdir '/' exp '/fmri/bined_pattern_regression/' timeUnit '/roi/' froidir '/coulpingF_sl_bin' num2str(binSize) 'lag' num2str(min(lags)) '-' num2str(max(lags)) '_fwemax.nii']);

    save([expdir '/' exp '/fmri/bined_pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_bin' num2str(binSize) 'lag'  lag'  num2str(min(lags)) '-' num2str(max(lags)) '_stats' ],'couple_sig_fdr','couple_sig_fwe','couple_sig_fwemax','F','p','lags','rnames');
    
    clear  couple_sig_fwe  couple_sig_fdr  couple_sig_fwemax
end
