clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

% -4 for merlin, -3 for sherlock
lags=-9:-4;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

for ei=1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LLminus1L_lag0-0'],'couplingz','rnames','keptT','b');
    cp_ll=squeeze(nanmean(couplingz,4));

    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SLminus1L_lag' num2str(min(lags)) '-' num2str(max(lags))],'b','couplingz','rnames','keptT','lags');
    cp_sl=squeeze(couplingz);

    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SLfake_lag' num2str(min(lags)) '-' num2str(max(lags))],'couplingz');
    cp_slf=squeeze(couplingz);
  
       load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag-10-10_stats' ],'couple_sig_fdr');
    for sp=1:size(couplingz,3);
         herd(:,sp)=corr_col(cp_sl(:,keptT,sp)',cp_ll(:,keptT,sp)');
        herd_null(:,sp)=corr_col(cp_slf(:,keptT,sp)',cp_ll(:,keptT,sp)');
    end
    
    % save herd value
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames);
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, herd]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '.nii']);
    
    herdz=real(atanh(herd));
    herdz_m=mean(herdz,2);
    herdz_null=real(atanh(herd_null));
    for ri=1:size(cp_sl,1);
        [~,herd_p(ri,1)]=ttest(herdz(ri,:)'-herdz_null(ri,:)',0,'Tail','right');
    end
 
    % test herd effect within regions showing significant coupling
    herd_sig_fdr=zeros(size(rnames));
    herd_sig_fdr(couple_sig_fdr==1)=fdr0(herd_p(couple_sig_fdr==1),0.05);
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_'  exp '_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'rnames','lags','herd','herd_p','herd_sig_fdr','herd_null');
rnames(herd_sig_fdr==1)

    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(herd_sig_fdr(:,1)==1));
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, herd(herd_sig_fdr(:,1)==1,1)]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr.nii']);
    
    clear herd herd_null herd_p herd_sig_fdr
end
