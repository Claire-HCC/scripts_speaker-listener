% flirt -in "aal_ROI_MNI_V4_2mm.nii" -ref "/scratch/claire/speaker-listener/roi_mask/MNI152NLin2009cAsym_3x3x4mm_brain.nii" -out "aal_ROI_MNI_V4_3x3x4mm.nii"  -applyxfm -usesqform -interp nearestneighbour
% fslchfiletype NIFTI aal_ROI_MNI_V4_3x3x4mm.nii aal_ROI_MNI_V4_3x3x4mm
clear all
loc='mypc';
set_parameters;

for ei=1:11;
    exp=exp_parameters.experiments{ei};
    nii=load_nii([expdir exp '\fmri\temporal\circularlagcorr\tr\vox\LL_leave1out\isc_r.nii']);
    
    roimask=load_nii([expdir '/roi_mask/aal/nii/Angular_L.nii']);
    roimask=roimask.img;
    
    nii.img(roimask==0)=0;
    mx=max(nii.img(:));
    
    [xv, yv, zv]=ind2sub(size(nii.img),find(nii.img==mx));
    cor=cor2mni([xv yv zv],[nii.hdr.hist.srow_x; nii.hdr.hist.srow_y ;nii.hdr.hist.srow_z; 0 0 0 1]);
    exp_parameters.isc_peak_Angular_L{ei}=cor;
end