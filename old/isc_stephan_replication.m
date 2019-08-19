clear all
close all
set_parameters

lag=20;

role='listener';
tic

mask=  load_nii([expdir 'roi_mask/MNI152_T1_3mm_brain_mask.nii']);
mask_mni=mask.img(:);

for ei=1%:2;
    exp= experiments{ei};
    load([expdir exp '/fmri/mat/wholeBrain/'  exp '_speaker_subj01.mat']);
    epi_s=zscore(data')';
    
    subjects=cellstr(ls([expdir exp '/fmri/mat/wholeBrain/' exp '_listener*mat']));
    %  fsize=[10 30];
    %  h=figure('position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize,'unit','centimeter');
    
    isc_peakR=zeros(size(epi_s,1),1,length(subjects));
    isc_peakT=zeros(size(epi_s,1),1,length(subjects));
    
    for si=1:length(subjects);
        subj=subjects{si};
        load([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj]);
        epi_l=zscore(data')';
        
        mask=find(sum(epi_s')'~=0 & sum(epi_l')'~=0 & mask_mni~=0);
        
        for vi = 1:length(mask);
            v=mask(vi);
            
            [r,lags]=xcorr(epi_l(v,:),epi_s(v,:),lag,'coeff'); %!!!
            pk=max(r);
            pt=lags(r==pk); pt=pt(1);
            isc_peakR(v,1,si)=pk;
            isc_peakT(v,1,si)=pt;
        end
        
        nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
        data=isc_peakR(:,1,si);
        data((end+1):(61*73*61),:)=0;
        nii.img=reshape(data,61,73,61);
        save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/isc_' subj '_peakR.nii']);
        
        nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
        data=isc_peakT(:,1,si);
        data((end+1):(61*73*61),:)=0;
        nii.img=reshape(data,61,73,61);
        save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/isc_' subj '_peakT.nii']);
        
    end
end

save([expdir exp '/fmri/mat/wholeBrain/isc_' role '_peakR.mat' ],'isc_peakR');
save([expdir exp '/fmri/mat/wholeBrain/isc_' role '_peakT.mat' ],'isc_peakT');


nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
data=mean(isc_peakR,3);
data((end+1):(61*73*61),:)=0;
nii.img=reshape(data,61,73,61);
save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/isc_' role '_peakR.nii']);

nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
data=mean(isc_peakT,3);
data((end+1):(61*73*61),:)=0;
nii.img=reshape(data,61,73,61);
save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/isc_' role '_peakT.nii']);

thr=0.2;
nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
data=sum((isc_peakR>0.2),3);
data((end+1):(61*73*61),:)=0;
nii.img=reshape(data,61,73,61);
save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/isc_' role '_aboveThrSubjn.nii']);