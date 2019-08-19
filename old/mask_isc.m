clear all
close all
set_parameters



for ei=1%:2;
    exp= experiments{ei};
    
    load([ expdir exp '/fmri/mat/wholeBrain/isc_listener_peakR.mat']);
    isc=mean(isc_peakR,3);
    
    load([expdir '/roi_mask/mat/' exp '_subjsMask.mat']);
    mask_subjs=mask;
    
    mask=zeros(voxn,1);
    
    i=intersect(find(mask_subjs~=0), find(isc >=0.12));
    mask(i)=1;
    
    save([expdir '/roi_mask/mat/isc_thr12_' exp '.mat'],'mask');
    
    nii=mat2nii(mask);
    save_nii(nii,[expdir '/roi_mask/nifti/isc_thr12_' exp '.nii']);
end
