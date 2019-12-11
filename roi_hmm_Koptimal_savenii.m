clear all
close all
set_parameters
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

exp='merlin';
tableK=readtable([expdir   exp '/fmri/hmm/Ks_optimal.csv'])
rnames=table2array(tableK(:,1));
Ks=table2array(tableK(:,2));

roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames,'UniformOutput',0);
ris=(cellfun(@isempty,roi_table_inds)==0);

Ks=Ks(ris);
roi_ids=cell2mat(roi_table.id(cell2mat(roi_table_inds(ris))));
nii=roiTable2wholeBrainNii_mor([roi_ids, Ks]);
save_nii(nii,[expdir   exp '/fmri/hmm/Ks_optimal.nii']);
