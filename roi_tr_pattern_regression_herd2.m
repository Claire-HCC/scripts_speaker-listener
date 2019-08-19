clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

lags=-4;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

for ei=3;%1:4;
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
    
    herd(:,1)=corr_col(cp_sl(:,keptT)',cp_ll(:,keptT)');
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_'  exp '_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'b_sl','b_ll','rnames','lags','herd');
    
    % save herd value
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames);
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, herd]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '.nii']);
    
    % stats using phase_andomized sl
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm'],'couplingz_perm');
    iters=size(couplingz_perm,3);
    cp_sl_perm=couplingz_perm;
    for iter=1:iters;
        herd_perm(:,iter)=corr_col(cp_sl_perm(:,keptT,iter)',cp_ll(:,keptT)');
    end
    herd_p=sum(herd<herd_perm,2)/iters;
 
    % stats only within regions showing significant sl coupling
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_stats' ],'couple_sig_fdr','couple_sig_fwe','couple_sig_fwemax');
    couple_sig=couple_sig_fwemax;
    herd_sig_fdr=zeros(length(rnames),1);   
    herd_sig_fdr(couple_sig==1,1)=fdr0(herd_p(couple_sig==1,1),0.05);
    rnames(herd_sig_fdr(:,1)==1)
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_'  exp '_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'b_sl','b_ll','herd_sig_fdr','rnames','lags','herd','herd_p','herd_perm','cp_sl','cp_ll','keptT');
    
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(herd_sig_fdr(:,1)==1));
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, herd(herd_sig_fdr(:,1)==1,1)]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr.nii']);
    
    clear herd herd_perm
end