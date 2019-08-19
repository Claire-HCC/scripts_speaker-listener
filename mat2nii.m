function nii=mat2nii(data);
set_parameters;
nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3x3x4mm_brain.nii']);
data=mean(data,3);
data((end+1):(voxn),:)=NaN;
nii.img=reshape(data,volsize(1),volsize(2),volsize(3));
end