clear all
close all
loc='mypc';
set_parameters

rnames={'HG_L','PMC_L','vPCUN','vmPFC','aANG_L','STC_L','pIFG_L'};
load([ expdir '/roi_mask/' 'merlin_groupAvtive_mask.mat'],'roimask');

role1='speaker01';
role2='listenerZscoreMean';
relation='SL';
lags=[-10:10];
p_thr=0.05;

for ei=1%:2;
    exp=experiments{ei};
    %  subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/' role '*.mat']));
    % subjects= subjects(~cellfun(@isempty,(regexp(subjects,'[0-9]*'))));
    % subjects=cellfun(@(x) strrep(x,'.mat',''),subjects,'UniformOutput',0);
    
    for ri=3%1:length(rnames);
        
        rname=rnames{ri};
        
        %     for si=1%:length(subjects);
        
        %   subj=subjects{si};
        
        load([expdir experiments{ei} '/fmri/mat/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_'  rname '_perm.mat'],'keptvox','peakR_perm');
        peakR_perm=squeeze(peakR_perm);
        
        peakR=load_nii([expdir experiments{ei} '/fmri/nii/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_'  num2str(min(lags)) '-' num2str(max(lags)) '_' rname '_peakR.nii']);
        peakT=load_nii([expdir experiments{ei} '/fmri/nii/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_'  num2str(min(lags)) '-' num2str(max(lags)) '_' rname  '_peakT.nii']);
        r=peakR.img(:);
        r=r(keptvox);
        t=peakT.img(:);
        t=t(keptvox);
        
        
        p=  bsxfun(@gt, peakR_perm', r');
        p=sum(p)/size(peakR_perm,2);
        p_3d=nan(voxn,1);
        p_3d(keptvox)=p;

        p_FDR_mask=fdr0(p,p_thr);
        p_FDR_mask_3d=nan(voxn,1);
        p_FDR_mask_3d(keptvox)=p_FDR_mask;
        
        
        p_FWE_mask=p<(p_thr/(size(peakR_perm,1)));
        p_FWE_mask_3d=nan(voxn,1);
        p_FWE_mask_3d(keptvox)=p_FWE_mask;
        
        
        p_3d_nii=peakR;
        p_3d_nii.img=p_3d;
        save_nii(p_3d_nii,[expdir exp '/fmri/nii/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_'  num2str(min(lags)) '-' num2str(max(lags))  '_' rname '_p.nii']);
        
        r_FWE=peakR;
        r_FWE.img(p_FWE_mask_3d==0)=NaN;
        t_FWE=peakT;
        t_FWE.img(p_FWE_mask_3d==0)=NaN;
        
        save_nii(t_FWE,[expdir exp '/fmri/nii/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_'  num2str(min(lags)) '-' num2str(max(lags))  '_' rname '_peakT_p' num2str(p_thr) 'FWE.nii']);
        save_nii(r_FWE,[expdir exp '/fmri/nii/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_'  num2str(min(lags)) '-' num2str(max(lags))  '_' rname '_peakR_p' num2str(p_thr) 'FWE.nii']);
        
        r_FDR=peakR;
        r_FDR.img(p_FDR_mask_3d==0)=NaN;
        t_FDR=peakT;
        t_FDR.img(p_FDR_mask_3d==0)=NaN;
        save_nii(t_FDR,[expdir exp '/fmri/nii/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_'  num2str(min(lags)) '-' num2str(max(lags))  '_' rname '_peakT_p' num2str(p_thr) 'FDR.nii']);
        save_nii(r_FDR,[expdir exp '/fmri/nii/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_'  num2str(min(lags)) '-' num2str(max(lags))  '_' rname '_peakR_p' num2str(p_thr) 'FDR.nii']);
        
    end
    
    % clear peakR_perm peakT peakR
    %    end
end





