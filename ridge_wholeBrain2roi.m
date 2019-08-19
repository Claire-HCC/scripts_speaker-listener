clear all
close all
set_parameters


role='listener';
lags=[-10:10]; %scan

rtable=readtable([expdir 'roi_mask/roi_id_region.txt'],'Delimiter',' ');
rois_selected={'vPCUN'};
rids=rtable.id(find(contains(rtable.region,rois_selected)));

mask=  load_nii([expdir 'roi_mask/MNI152_T1_3mm_brain_mask.nii']);
mask_mni=mask.img(:);

for ei=1%:2;
    exp=experiments{ei};
    
    load([expdir experiments{ei} '/fmri/mat/wholeBrain/speaker01.mat']);
    speaker=data;
    
    for ri=1%:length(rids);
        rid=rids{ri};
        rname=rois_selected{ri};
        load([expdir '/roi_mask/mat/' rid '.mat']);
        rmask=data;
        
        speaker_roi=mean(speaker(rmask,:));
        speaker_roi=zscore(speaker_roi')';
        
        %  subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/*' role '*subj*mat']));
        %  subjects= subjects(~cellfun(@isempty,(regexp(subjects,'[0-9]*'))));
        subjects={'listenerZscoreMean'};
        
         B_slm=zeros(voxn,length(lags));
         F_slm=zeros(voxn,1);
         p_slm=zeros(voxn,1);
        
        for si=1:length(subjects);
            subj=subjects{si};
            subjID=char(regexp(subj,'subj..','match'));
            
            load([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj]);
            listener=zscore(data')';
            
            
            vox_selected=find(sum(listener')'~=0 & mask_mni~=0);
            speaker_roi_mat=repmat(speaker_roi,length( vox_selected),1);
            
            % use future speaker data to explain current listener data
            [B,F,p] = lag_ridgeregression3(listener(vox_selected,:), speaker_roi_mat, [lags(1) lags(end)]);
            
        end
        
        
        B_slm(vox_selected,:)=B;
        F_slm(vox_selected)=F;
        p_slm(vox_selected)=p;
        
        
        clear B F p
        
        B=B_slm; F=F_slm; p=p_slm;
        save([expdir experiments{ei} '/fmri/mat/wholeBrain/isfc_seed/isfc_seed_slm_' rname '.mat'],'B','F','p');
    end
end



nii=mat2nii(F);
save_nii(nii,[expdir experiments{ei} '/fmri/nii/wholeBrain/isfc_seed/isfc_seed_slm_' rname '_F.nii']);

p_temp=sort(p(p~=0));
FDR = mafdr(p_temp);
fdri=max(find(FDR<=0.01));
p_FDR=p_temp(fdri);
F_thr=finv(1-p_FDR,21,585-22)


