



load('Y:\claire\ABC\fMRI\timeseries\tr\wholeBrain\groupAll.mat');

for s=2:25;
    mat=zeros(voxn,size(gdata,2));
    mat(keptvox,:)=gdata(:,:,s);
    
    nii=mat2nii(mat);
    save_nii(nii,sprintf('Y:/claire/speaker-listener/ABC/fmri/timeseries/tr_uncropped/wholeBrain/old/listener%02d.nii',s))
end
