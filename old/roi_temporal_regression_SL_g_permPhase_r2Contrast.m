clear all
close all

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
roi_ids=cell2mat(roi_table.id);
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
lags1=[-10:-4];
lags2=[4:10];

for ei=4;%1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags1)) '-' num2str(max(lags1)) '_permPhase.mat' ],'b','r2','rnames');
    r2_perm1=r2;
    b_perm1=b;
    
    load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags2)) '-' num2str(max(lags2)) '_permPhase.mat' ],'b','r2','rnames');
    r2_perm2=r2;
    b_perm2=b;
    
    r2_perm=r2_perm2-r2_perm1;
    
    load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags1)) '-' num2str(max(lags1)) '.mat' ],'b','r2','rnames');
    r2_1=r2;
    b_1=b;
    
    load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags2)) '-' num2str(max(lags2)) '.mat' ],'b','r2','rnames');
    r2_2=r2;
    b_2=b;
   r2= r2_2-r2_1;
    
    p=min(mean(r2_perm>r2,2) , mean(r2_perm<r2,2));
    p(isnan(r2))=NaN;
    p=squeeze(p);
    pfdr=nan(size(p));
    [~,~,pfdr(~isnan(p(:)))]=fdr(p(~isnan(p)));
    pfdr=reshape(pfdr,size(p));
    
    pfwe=p*sum(~isnan(p(:)));
    
    rnames(pfdr<.025)
    
    %         save([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_r2Contrast.mat' ],'r2','sig_fdr','p','pfdr','pfwe','b');
    %
    %         rnames(pfdr<.05)
    %         nii=roiTable2wholeBrainNii_mor([roi_ids(pfdr<.05),  r2(pfdr<.05)]);
    %         save_nii(nii,[expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_g/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_r2Contrast_pfdr.nii']);
    %
    %  %      rnames(pfwe<.05)
    %         nii=roiTable2wholeBrainNii_mor([roi_ids( pfwe<.05),  r2( pfwe<.05)]);
    %         save_nii(nii,[expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_g/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_r2Contrast_ pfwe.nii']);
end


