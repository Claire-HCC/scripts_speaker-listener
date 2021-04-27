
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
    
    if ei==1;
        nii=load_nii([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/coulpingPeakLag_sl_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr.nii']);
        nii.img=(~isnan(nii.img));
    else
        nii2=load_nii([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/coulpingPeakLag_sl_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr.nii']);
        nii.img=nii.img+(~isnan(nii2.img));
    end 
end

save_nii(nii,[expdir '/roi_mask/coulpingcoulpingPeakLag_sl_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr_overlap.nii']);
nii.img=(nii.img==4);
save_nii(nii,[expdir '/roi_mask/coulpingcoulpingPeakLag_sl_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr_overlap4.nii']);
    
