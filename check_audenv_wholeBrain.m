clear all
close all
set_parameters

mask=  load_nii([expdir 'roi_mask/MNI152_T1_3mm_brain_mask.nii']);
mask_mni=mask.img(:);

role='speaker';
lag=6;

for ei=1:2;
    exp=experiments{ei};
    
    load([expdir experiments{ei} '/sound/' experiments{ei} '_audenv.mat' ]);
    aud=zscore(aud);
    
    subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/*' role '*subj*mat']));
    
    for si=1:length(subjects);
        subj=subjects{si};
        load([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj]);
        epi=zscore(data')';
        
        mask=find(sum(data')'~=0 & mask_mni~=0);
        for vi = 1:length(mask);
            v=mask(vi);
            %  epi(v,:)=despike(epi(v,:)',[-3 3],5);
            % ti=find(epi(v,:) > -3 & epi(v,:) <3);
            
            [r,lags]=xcorr(epi(v,:),aud(:),lag,'coeff');
            pk=max(r);
            pt=lags(r==pk); pt=pt(1);
            audenv_peakR(v,1,si)=pk;
            audenv_peakT(v,1,si)=pt;
        end
        
        
        nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
        data=audenv_peakR(:,1,si);
        data((end+1):(61*73*61),:)=0;
        nii.img=reshape(data,61,73,61);
        save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/audenv_' subj '_peakR.nii']);
        
        nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
        data=audenv_peakT(:,1,si);
        data((end+1):(61*73*61),:)=0;
        nii.img=reshape(data,61,73,61);
        save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/audenv_' subj '_peakT.nii']);
        
        
    end
    save([expdir exp '/fmri/mat/wholeBrain/' exp '_audenv_' role '_peakR.mat' ],'audenv_peakR');
    save([expdir exp '/fmri/mat/wholeBrain/' exp '_audenv_' role '_peakT.mat' ],'audenv_peakT');
    
    
    nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
    data=mean(audenv_peakR,3);
    data((end+1):(61*73*61),:)=0;
    nii.img=reshape(data,61,73,61);
    save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/audenv_' role '_peakR.nii']);
    
    nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
    data=mean(audenv_peakT,3);
    data((end+1):(61*73*61),:)=0;
    nii.img=reshape(data,61,73,61);
    save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/audenv_' role '_peakT.nii']);
    
end


