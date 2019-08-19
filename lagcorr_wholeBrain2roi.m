
% role1 roi to role2 wholebrain
% negative peak T means role1 precedes, positive peakT means role2 precedes
clear all
close all
loc='mypc';
set_parameters

role1='listener01';
role2='listener01_others';
relation='LL';
lags=[-10:10]; %scan

rnames={'HG_L','PMC_L','vPCUN','vmPFC','aANG_L','STC_L','pIFG_L'};

mask=  load_nii([expdir 'roi_mask/rspm8_gray0.5_merlin_groupAvtive_mask.nii']);
mask=mask.img(:);

for ei=1%:2;%1:2;
    exp=experiments{ei};
    
    load([expdir experiments{ei} '/fmri/mat/wholeBrain/' role1 '.mat']);
    data1=data;
    
    for ri=3%1:length(rnames);
        rname=rnames{ri};
        load([expdir '/roi_mask/mat/' rname '.mat']);
        
        role1_roi=mean(data1(roimask==1,:));
        role1_roi=zscore(role1_roi,0,2);
        
        %  subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/*' role '*.mat']));
        % subjects= subjects(~cellfun(@isempty,(regexp(subjects,'[0-9]*'))))
        
        
        % compute listener-speaker isfc for each listener
        %   for si = 1%:size(listeners,3);
        % choose listener and take zscore of timeseries in each voxel.
        %   [filepath,subj,ext] = fileparts(subjects{si});
        load([expdir exp '/fmri/mat/wholeBrain/' role2 '.mat']);
        data2=data;
        
        keptvox=find(mask==1);
        
        role1_roi_mat=repmat(role1_roi,length(keptvox),1);
        
        i=400:500;
        for si=1:size(data,3);
            lagcc_temp= lagcorr_claire(role1_roi_mat(:,i)',data(keptvox,i,si)',lags)';
            [peakR(:,1,si) peakT(:,1,si)]=max((lagcc_temp)');
            peakT(:,1,si)=lags(peakT(:,1,si))';
        end
        
        save([expdir experiments{ei} '/fmri/mat/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_' num2str(min(lags)) '-' num2str(max(lags)) '_'  rname '.mat'],'peakT','peakR','keptvox');
        
        
        peakR_nii=nan(voxn,1);
        peakR_nii(keptvox,:)=peakR;
        nii=mat2nii(peakR_nii);
        save_nii(nii,[expdir experiments{ei} '/fmri/nii/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_' num2str(min(lags)) '-' num2str(max(lags)) '_'  rname '_peakR.nii']);
        
        peakT_nii=nan(voxn,1);
        peakT_nii(keptvox,:)=peakT;
        nii=mat2nii(peakT_nii);
        save_nii(nii,[expdir experiments{ei} '/fmri/nii/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_' num2str(min(lags)) '-' num2str(max(lags)) '_'  rname '_peakT.nii']);
        
        r_thr=0.2;
        peakT_nii=nan(voxn,1);
        peakT_nii(keptvox,:)=peakT;
        peakT_nii(peakR_nii<r_thr)=NaN;
        nii=mat2nii(peakT_nii);
        save_nii(nii,[expdir experiments{ei} '/fmri/nii/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_' num2str(min(lags)) '-' num2str(max(lags)) '_'  rname '_peakT_rthr' num2str(r_thr) '.nii']);
        
        %      end
        
    end
end

beep



