
clear all
close all

set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');

lags=-30:-4;
for ei=3;%1:4;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_stats.mat'],'r2_sig_fwe','rnames');
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lagSelection_stats'  ],'lags_peak');
    
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(r2_sig_fwe==1),'UniformOutput',0);
    roi_ids=cell2mat(roi_table.id(cell2mat(roi_table_inds)));
    nii=roiTable2wholeBrainNii_mor([roi_ids, lags_peak(r2_sig_fwe==1)]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lagSelectionPeak.nii']);
    
end
