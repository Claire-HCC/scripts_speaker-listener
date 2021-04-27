% the data is from Princeton Dataspace https://dataspace.princeton.edu/jspui/handle/88435/dsp015d86p269k.
clear all
close all
set_parameters

tr=1.5;
voln=272;
voxn=61*73*61;
% 12 s of neutral music and 3 s of silence preceded and 15 s of silence followed each playback in all conditions.

subjs_intact=1:36;
subjs_word=[1 3:4 10:14 16 20 23:28 30:32 34:35 37:51];

for si=1:length(subjs_intact);
    subj=subjs_intact(si);
    fname=ls([expdir '\fMRI\roi_mask\pieman_TRW\sub-' num2str(si) '\func\sub-' num2str(subj) '-task-intact*.nii']);;
    nii=load_nii([expdir '\roi_mask\pieman_TRW\sub-' num2str(si) '\func\' fname]);
    nii.img=nii.img(:,:,:,(crop_startn+1):(crop_startn+voln));
    gdata(:,:,si)=reshape(nii.img,voxn,voln);
end
keptvox=find((sum(sum(gdata==0,2),3)==0));
gdata=gdata(keptvox,:,:);
save([expdir '\roi_mask\pieman_TRW\allSubjs_intact.mat'],'gdata','keptvox','-v7.3');

clear gdata
for si=1:length(subjs_word);
    subj=subjs_word(si);
    fname=ls([expdir '\roi_mask\pieman_TRW\sub-' num2str(subj) '\func\sub-' num2str(subj) '-task-word*.nii']);;
    nii=load_nii([expdir '\roi_mask\pieman_TRW\sub-' num2str(subj) '\func\' fname]);
    nii.img=nii.img(:,:,:,(crop_startn+1):(crop_startn+voln));
    gdata(:,:,si)=reshape(nii.img,voxn,voln);
end
keptvox=find((sum(sum(gdata==0,2),3)==0));
gdata=gdata(keptvox,:,:);
save([expdir '\roi_mask\pieman_TRW\allSubjs_word.mat'],'gdata','keptvox','-v7.3');
