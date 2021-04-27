clear all
close all
set_parameters


role='listener';
lag=20;

rtable=readtable([expdir 'roi_mask/roi_id_region.txt'],'Delimiter',' ');
rois_selected={'HG_L','Precentral_L','PCC'};
rids=rtable.id(find(contains(rtable.region,rois_selected)));

mask=  load_nii([expdir 'roi_mask/MNI152_T1_3mm_brain_mask.nii']);
mask_mni=mask.img(:);

for ei=1%:2;
    exp=experiments{ei};
    
    load([expdir experiments{ei} '/fmri/mat/wholeBrain/' exp '_speaker_subj01.mat']);
    epi_s=zscore(data')';
    
    for ri=1:length(rids);
        rid=rids{ri};
        load([expdir '/roi_mask/mat/' rid '.mat']);
        rmask=data;
        
        epi_s_masked=mean(epi_s(rmask,:));
        
        subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/*' role '*subj*mat']));
        
        for si=1:length(subjects);
            subj=subjects{si};
            subjID=char(regexp(subj,'subj..','match'));
            load([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj]);
            epi_l=zscore(data')';
            
            mask=find(sum(epi_l')'~=0 & mask_mni~=0);
            
            for vi = 1:length(mask);
                v=mask(vi);
                
                  [r,lags,t] = migram(double(epi_l(v,:)),double(epi_s_masked),20,length(epi_s_masked),0,'norm');
                pk=max(r);
                pt=lags(r==pk); pt=pt(1);
                roi2wholeBrain_peakR(v,1,si)=pk;
                roi2wholeBrain_peakT(v,1,si)=pt;
            end
            
            
            nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
            data=roi2wholeBrain_peakR(:,1,si);
            data((end+1):(61*73*61),:)=0;
            nii.img=reshape(data,61,73,61);
            save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/roi2wholeBrain_mi_' rois_selected{ri} '_' subjID   '_peakR.nii']);
            
            nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
            data=roi2wholeBrain_peakT(:,1,si);
            data((end+1):(61*73*61),:)=0;
            nii.img=reshape(data,61,73,61);
            save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/roi2wholeBrain_mi_' rois_selected{ri} '_' subjID '_peakT.nii']);

        end
        save([expdir exp '/fmri/mat/wholeBrain/roi2wholeBrain_mi_' rois_selected{ri}   '_peakR.mat' ],'roi2wholeBrain_peakR');
        save([expdir exp '/fmri/mat/wholeBrain/roi2wholeBrain_mi_' rois_selected{ri}   '_peakT.mat' ],'roi2wholeBrain_peakT');
        
        nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
        data=mean(roi2wholeBrain_peakR,3);
        data((end+1):(61*73*61),:)=0;
        nii.img=reshape(data,61,73,61);
        save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/roi2wholeBrain_mi_' rois_selected{ri}   '_peakR.nii']);
        
        nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
        data=mean(roi2wholeBrain_peakT,3);
        data((end+1):(61*73*61),:)=0;
        nii.img=reshape(data,61,73,61);
        save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/roi2wholeBrain_mi_' rois_selected{ri}  '_peakT.nii']);
        
        thr=0.03;
        nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
        data=sum((roi2wholeBrain_peakR>0.2),3);
        data((end+1):(61*73*61),:)=0;
        nii.img=reshape(data,61,73,61);
        save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/roi2wholeBrain_mi_' rois_selected{ri}  '_aboveThrSubjn.nii']);
        
    end
end


