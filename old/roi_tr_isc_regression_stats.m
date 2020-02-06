clear all;
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
lags=-4:4;

for ei=3;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F_sl','b_sl','rnames_sl');
    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],'F_sl_perm','b_sl_perm');
    
    %% significant coupling
    p=sum(F_sl<F_sl_perm,3)/1000;
    couple_sig_fwe(:,1)=(p<(0.05/length(rnames_sl)));
    couple_sig_fdr(:,1)=fdr0(p,0.05);
    couple_sig_fwemax=(F_sl>quantile(max(squeeze(F_sl_perm)),0.95));
    
    %     % save F_sl
    %     roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames_sl);
    %     roi_ids=cell2mat(roi_table.id(roi_table_inds));
    %     nii=roiTable2wholeBrainNii_mor([roi_ids, F_sl]);
    %     save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/coulpingF_sl_' num2str(min(lags)) '-' num2str(max(lags)) '.nii']);
    %
    %     roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames_sl(couple_sig_fwe(:,1)==1));
    %     roi_ids=cell2mat(roi_table.id(roi_table_inds));
    %     nii=roiTable2wholeBrainNii_mor([roi_ids, F_sl(couple_sig_fwe(:,1)==1,1)]);
    %     save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/coulpingF_sl_' num2str(min(lags)) '-' num2str(max(lags)) '_fwe.nii']);
    %
    %     roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames_sl(couple_sig_fdr(:,1)==1));
    %     roi_ids=cell2mat(roi_table.id(roi_table_inds));
    %     nii=roiTable2wholeBrainNii_mor([roi_ids, F_sl(couple_sig_fdr(:,1)==1,1)]);
    %     save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/coulpingF_sl_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr.nii']);
    %
    %     roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames_sl(couple_sig_fwemax(:,1)==1));
    %     roi_ids=cell2mat(roi_table.id(roi_table_inds));
    %     nii=roiTable2wholeBrainNii_mor([roi_ids, F_sl(couple_sig_fwemax(:,1)==1,1)]);
    %     save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/coulpingF_sl_' num2str(min(lags)) '-' num2str(max(lags)) '_fwemax.nii']);
    
    %% peak lags
    % peak lag in each region
    
    
    % compute p-value (the first beta value corresponds to intercept)
    b_p=sum(b_sl_perm(:,2:end,:)>b_sl(:,2:end),3)/size(b_sl_perm,3);
    sig_b=b_p<(0.05/size(b_sl,1)*size(b_sl,2));
    
    % focus on beta values that are significant
    b_temp=b_sl(:,2:end);
    b_temp(sig_b~=1)=0;
    
    % find peak lags
    [v col]=max(b_temp(:,2:end),[],2);
    col(sum(sig_b,2)==0)=nan;
    bpt=nan(size(col));
    bpt(~isnan(col),1)=lags(col(~isnan(col)));
    
    %     figure; hist(bpt(couple_sig_fdr(:,1)==1),21); title(exp);
    %     xlim([-10 10]);
    %     ylim([0 12])
    
    rnames_sl(~isnan(bpt) & (couple_sig_fdr(:,1)==1))
    % save peak b_sl lag
    % to display this in xjview
    %  colormap([gray(64);[flipud(hot(32)); winter(32)]])
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames_sl);
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, bpt]);
    save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/coulpingPeakLag_sl_' num2str(min(lags)) '-' num2str(max(lags)) '.nii']);
    
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames_sl(couple_sig_fdr(:,1)==1));
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, bpt(couple_sig_fdr(:,1)==1,1)]);
    save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/coulpingPeakLag_sl_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr.nii']);
    
    save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_stats' ],'couple_sig_fdr','couple_sig_fwe','couple_sig_fwemax','F_sl','p','lags','rnames_sl','bpt');
    
    clear  couple_sig_fwe  couple_sig_fdr  couple_sig_fwemax bpt
end
