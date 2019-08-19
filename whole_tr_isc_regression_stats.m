clear all;
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
lags=-10:10;

for ei=3;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F_sl','rnames_sl','lags_sl');
    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],'F_sl_perm');
    
    p_sl=sum(F_sl<F_sl_perm,3)/1000;
    couple_sl_sig_fwe=(p_sl<(0.05/length(F_sl)));
    couple_sl_sig_fdr=fdr0(p_sl,0.05);
    
    save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_stats' ],'couple_sl_sig_fdr','couple_sl_sig_fwe','p_sl','lags','keptvox');
    
    
    
    %% peak lags
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
    
    
    
    %     rnames_sl(~isnan(bpt) & (couple_sig_fdr(:,1)==1))
    %     % save peak b_sl lag
    %     % to display this in xjview
    %     %  colormap([gray(64);[flipud(hot(32)); winter(32)]])
    %     roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames_sl);
    %     roi_ids=cell2mat(roi_table.id(roi_table_inds));
    %     nii=roiTable2wholeBrainNii_mor([roi_ids, bpt]);
    %     save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' froidir '/coulpingPeakLag_sl_' num2str(min(lags)) '-' num2str(max(lags)) '.nii']);
    %
    %     roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames_sl(couple_sig_fdr(:,1)==1));
    %     roi_ids=cell2mat(roi_table.id(roi_table_inds));
    %     nii=roiTable2wholeBrainNii_mor([roi_ids, bpt(couple_sig_fdr(:,1)==1,1)]);
    %     save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' froidir '/coulpingPeakLag_sl_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr.nii']);
    
    save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' froidir '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_stats' ],'couple_sig_fdr','couple_sig_fwe','F_sl','p','lags','rnames_sl','bpt');
    
    clear  couple_sig_fwe  couple_sig_fdr  couple_sig_fwemax bpt
end
