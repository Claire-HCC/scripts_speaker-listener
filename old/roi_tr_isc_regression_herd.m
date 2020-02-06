
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

lags=-4:4;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');


for ei=3%:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags))],'rnames','keptT_ll','b_ll','r_ll','Y_ll');
    b_ll=b_ll(:,2:end,:);
    
    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags))],'b_sl','rnames_sl','keptT_sl','lags_sl','r_sl','Y_sl');
    b_sl=b_sl(:,2:end);
    rnames=rnames_sl;
    ketpT=keptT_sl;
    lags=lags_sl;
    
    %    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) 'stats'],'coulping_sl_sig_fdr','herd_sig_fwe');
    
    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm'],'r_sl_perm','Y_sl_perm');
    iters=size(r_sl_perm,3);
    
    %  percent absolute error
    y_sl=Y_sl+r_sl;
    pabsr_sl=(abs(r_sl));
    
    y_ll=Y_ll+r_ll;
    pabsr_ll=nanmean(abs(r_ll),3);
    
    herd(:,1)=corr_col(pabsr_sl',pabsr_ll');
    
    y_sl_perm=Y_sl_perm+r_sl_perm;
    pabsr_sl_perm=(abs(r_sl_perm));
    for iter=1:iters;
        herd_perm(:,iter)=corr_col(pabsr_sl_perm(:,:,iter)',pabsr_ll');
    end
    
    herd_p(:,1)=sum((herd_perm>repmat(herd,1,iters)),2)/iters;
    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_stats' ],'couple_sl_sig_fdr','couple_sl_sig_fwe','couple_sig_fwemax');
    couple_sig=couple_sl_sig_fwe;
    herd_sig_fdr=zeros(length(rnames),1);
    herd_sig_fdr(couple_sig==1,1)=fdr0(herd_p(couple_sig==1,1),0.05);
    rnames(herd_sig_fdr(:,1)==1)
    
    herd_sig_fwe=zeros(length(rnames),1);
    herd_sig_fwe(couple_sig==1,1)=(herd_p(couple_sig==1,1)/sum(couple_sig)<0.05);
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(herd_sig_fwe(:,1)==1));
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, herd(herd_sig_fwe(:,1)==1,1)]);
    save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '_fwe.nii']);
    
    save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/herd_'  exp '_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'b_sl','b_ll','rnames','lags','herd','herd_p','herd_perm','pabsr_sl','pabsr_ll','pabsr_sl_perm');
    
    clear herd herd_perm
end
