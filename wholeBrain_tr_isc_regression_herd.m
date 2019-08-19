clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
lags=-4:4;

for ei=3%1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags))],'keptT_ll','b_ll','r_ll','Y_ll');
    b_ll=b_ll(:,2:end,:);
    
    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags))],'b_sl','keptT_sl','lags_sl','r_sl','Y_sl');
    b_sl=b_sl(:,2:end);
    ketpT=keptT_sl;
    lags=lags_sl;
    
    %  percent absolute error
    y_sl=Y_sl+r_sl;
    pabsr_sl=(abs(r_sl./y_sl));
    
    y_ll=Y_ll+r_ll;
    pabsr_ll=nanmean(abs(r_ll./y_ll),3);
    
    herd(:,1)=corr_col(pabsr_sl',pabsr_ll');
    
    load([expdir '/roi_mask/group_fmri_mask_' exp '.mat'],'roimask');
    keptvox=find(roimask);
    mat=zeros(voxn,1);
    mat(keptvox)=herd;
    nii=mat2nii(mat);
    % save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/herd_lag' num2str(min(lags)) '-' num2str(max(lags)) '.nii']);

%     %    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) 'stats'],'coulping_sl_sig_fdr','herd_sig_fwe');
%     load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm'],'r_sl_perm','Y_sl_perm');
%     iters=size(r_sl_perm,3);
%     
%     y_sl_perm=Y_sl_perm+r_sl_perm;
%     pabsr_sl_perm=(abs(r_sl_perm./y_sl_perm));
%     for iter=1:iters;
%         herd_perm(:,iter)=corr_col(pabsr_sl_perm(:,:,iter)',pabsr_ll');
%     end
%     
%     herd_p(:,1)=sum((herd_perm>repmat(herd,1,iters)),2)/iters;
%     
%     herd_sig_fwe=(herd_p<(0.05/size(rnames_sl,1)));
%     
%     roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(herd_sig_fwe(:,1)==1));
%     roi_ids=cell2mat(roi_table.id(roi_table_inds));
%     nii=roiTable2wholeBrainNii_mor([roi_ids, herd(herd_sig_fwe(:,1)==1,1)]);
%     save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '_fwe.nii']);
%     
%     save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/herd_'  exp '_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'b_sl','b_ll','herd_sig_fwe','rnames','lags','herd','herd_p');
%     
end
