clear all
close all
set_parameters


mask=zeros()
for ei=1%:2;
    exp= experiments{ei};
    subjects=cellstr(ls([expdir exp '/fmri/mat/wholeBrain/' exp '_listener*mat']));
    
    
    mask=ones(voxn,1);
    
  for si=1:length(subjects);
        subj=subjects{si};
        load([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj]);
        epi=sum(data,1);
        i=find(epi==0 | isnan(epi));
        mask(i)=0;
  end
end

save([expdir '/roi_mask/mat/' exp '_subjsMask.mat'],'mask');

nii=mat2nii(mask);
save_nii(nii,[expdir '/roi_mask/nifti/' exp '_subjsMask.nii']);
