clear all;
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
lags=-10:10;

for ei=1:4;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F','b','rnames');
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],'F_perm','b_perm','rnames');
    
    %% significant coupling
    p=sum(F<F_perm,3)/1000;
    couple_sig_fwe(:,1)=(p<(0.05/length(rnames)));
    couple_sig_fdr(:,1)=fdr0(p,0.05);
    couple_sig_fwemax=(F>quantile(max(squeeze(F_perm)),0.95));
    
    %     % save F
    %     roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames);
    %     roi_ids=cell2mat(roi_table.id(roi_table_inds));
    %     nii=roiTable2wholeBrainNii_mor([roi_ids, F]);
    %     save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/coulpingF_sl_' num2str(min(lags)) '-' num2str(max(lags)) '.nii']);
    %
    %     roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(couple_sig_fwe(:,1)==1));
    %     roi_ids=cell2mat(roi_table.id(roi_table_inds));
    %     nii=roiTable2wholeBrainNii_mor([roi_ids, F(couple_sig_fwe(:,1)==1,1)]);
    %     save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/coulpingF_sl_' num2str(min(lags)) '-' num2str(max(lags)) '_fwe.nii']);
    %
    %     roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(couple_sig_fdr(:,1)==1));
    %     roi_ids=cell2mat(roi_table.id(roi_table_inds));
    %     nii=roiTable2wholeBrainNii_mor([roi_ids, F(couple_sig_fdr(:,1)==1,1)]);
    %     save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/coulpingF_sl_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr.nii']);
    %
    %     roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(couple_sig_fwemax(:,1)==1));
    %     roi_ids=cell2mat(roi_table.id(roi_table_inds));
    %     nii=roiTable2wholeBrainNii_mor([roi_ids, F(couple_sig_fwemax(:,1)==1,1)]);
    %     save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/coulpingF_sl_' num2str(min(lags)) '-' num2str(max(lags)) '_fwemax.nii']);
    
    %% peak lags
    % peak lag in each region
    
    
    % compute p-value (the first beta value corresponds to intercept)
    b_p=sum(b_perm(:,2:end,:)>b(:,2:end),3)/size(b_perm,3);
    sig_b=b_p<(0.05/size(b,1)*size(b,2));
    
    % focus on beta values that are significant
    b_temp=b(:,2:end);
    b_temp(sig_b~=1)=0;
    
    % find peak lags
    [v col]=max(b_temp(:,2:end),[],2);
    col(sum(sig_b,2)==0)=nan;
    bpt=nan(size(col));
    bpt(~isnan(col),1)=lags(col(~isnan(col)));
    
    %     figure; hist(bpt(couple_sig_fdr(:,1)==1),21); title(exp);
    %     xlim([-10 10]);
    %     ylim([0 12])
    
    rnames(~isnan(bpt) & (couple_sig_fdr(:,1)==1))
    % save peak b lag
    % to display this in xjview
    %  colormap([gray(64);[flipud(winter(32)); (hot(32)) ]])
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames);
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, bpt]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/coulpingPeakLag_sl_' num2str(min(lags)) '-' num2str(max(lags)) '.nii']);
    
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(couple_sig_fdr(:,1)==1));
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, bpt(couple_sig_fdr(:,1)==1,1)]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/coulpingPeakLag_sl_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr.nii']);
    
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_stats' ],'couple_sig_fdr','couple_sig_fwe','couple_sig_fwemax','F','p','lags','rnames','bpt');
    
    clear  couple_sig_fwe  couple_sig_fdr  couple_sig_fwemax bpt
end
