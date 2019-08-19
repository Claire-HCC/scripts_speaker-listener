function nii=mat2nii_3x3x3mm(data);
set_parameters;
volsize=[61 73 61];
voxn=prod(volsize);
nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
data=mean(data,3);
data((end+1):(voxn),:)=NaN;
nii.hdr.dime.dim(5)=size(data,2);
nii.hdr.dime.dim(1)=3+(nii.hdr.dime.dim(5)~=1);
nii.img=reshape(data,volsize(1),volsize(2),volsize(3),size(data,2));
end