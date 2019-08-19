% bined r2 vs std between subjects
clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

% -4 for merlin, -3 for sherlock
lags=-7:-4;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

for ei=1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_SL_minus1listener_lag' num2str(min(lags)) '-' num2str(max(lags))],'r2','rnames','keptT');
    r2_sl=r2;
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_SLfake_lag' num2str(min(lags)) '-' num2str(max(lags))],'r2');
    r2_slf=r2;
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_minus1listener_lag0-0'],'r2');
    r2_ll=nanmean(r2,4);
%         load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_lag0-0'],'r2');
%         r2_ll=nanmean(r2,3);
%         r2_ll=repmat(r2_ll,1,1,size(r2_sl,3));
%     
    % try dividing the story into early and late parts
    keptT=min(keptT):size(r2_sl,2);
    
    for sp=1:size(r2,3);
        for ri=1:size(r2,1)'
            herd(ri,sp)=corr(r2_sl(ri,keptT,sp)',r2_ll(ri,keptT,sp)','type','spearman');
            herd_null(ri,sp)=corr(r2_slf(ri,keptT,sp)',r2_ll(ri,keptT,sp)','type','spearman');
        end
    end
    
    herdz=real(atanh(herd));
    herdz_m=mean(herdz,2);
    herdz_null=real(atanh(herd_null));
    herdz_d=(herdz-herdz_null);
    herdz_d_m=nanmean(herdz_d')';
    [~,herd_p]=ttest([herdz-herdz_null]',0);
    
    herd_sig_fdr=zeros(size(herd_p));
    herd_sig_fdr=fdr0(herd_p,0.05);
    herd_sig_fdr_pos=zeros(size(herd_p));
    herd_sig_fdr_neg=zeros(size(herd_p));
    herd_sig_fdr_pos(herd_sig_fdr==1 & herdz_d_m'>0)=1;
    herd_sig_fdr_neg(herd_sig_fdr==1 & herdz_d_m'<0)=1;
    
    rnames(herd_sig_fdr_pos==1)
    % test herd effect within regions showing significant coupling
    % load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag-10-10_stats'],'couple_sig_fdr');
    % herd_sig_fdr(couple_sig_fdr==1)=fdr0(herd_p(couple_sig_fdr==1),0.05);
    
    herdm=mean(herd,2);
    herdm_null=nanmean(herd_null,2);
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames);
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, herdm]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '.nii']);
    
    herdm=mean(herd,2);
    herdm_null=nanmean(herd_null,2);
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames);
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, herdz_d_m]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herdzD_' num2str(min(lags)) '-' num2str(max(lags)) '.nii']);
    
    
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(herd_sig_fdr_pos==1));
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, herdz_d_m(herd_sig_fdr_pos==1)]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herdzD_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr_pos.nii']);
    
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(herd_sig_fdr_neg==1));
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, herdz_d_m(herd_sig_fdr_neg==1)]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herdzD_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr_neg.nii']);
    
    save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'lags','herd','herd_p','herdz_d','herd_sig_fdr_pos','herd_sig_fdr_neg','herd_null','rnames');
    
    
    clear herd herd_null herd_p herd_sig_fdr_pos herd_sig_fdr_neg herdz_d_m
end
