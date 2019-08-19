
clear all
close all
set_parameters

roles={'listenerZscoreMean','speaker'}; % 'speaker'

lags=-10:10;
p_thr=0.05;

for rolei=1%:length(roles);
    role=roles{rolei};
    
    for ei=1%:4;
        exp=experiments{ei};
        subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/' role '*.mat']));
        % subjects= subjects(~cellfun(@isempty,(regexp(subjects,'[0-9]*'))));
        subjects=cellfun(@(x) strrep(x,'.mat',''),subjects,'UniformOutput',0);
        
        for si=1:length(subjects);
            
            [filepath,subj,ext] = fileparts(subjects{si});
            
            load([expdir exp '/fmri/mat/wholeBrain/aud/audenv_' subj '_perm.mat']);
            peakR_perm=squeeze(peakR_perm);
            
            peakR=load_nii([expdir exp '/fmri/nii/wholeBrain/aud/audenv_' subj '_peakR.nii']);
            peakT=load_nii([expdir exp '/fmri/nii/wholeBrain/aud/audenv_' subj '_peakT.nii']);
            r=peakR.img(:);
            r=r(keptvox);
            t=peakT.img(:);
            t=t(keptvox);
            
            p=  bsxfun(@gt, peakR_perm', r');
            p=sum(p)/size(peakR_perm,2);
            p_3d=nan(voxn,1);
            p_3d(keptvox)=p;
            
            [p_sorted]=sort(p,'descend');
            [~,p_order]=ismember(p,p_sorted);
            p_FDR=p_order/(size(peakR_perm,1)+length(lags))*p_thr;
            p_FDR_mask=p<p_FDR;
            p_FDR_mask_3d=nan(voxn,1);
            p_FDR_mask_3d(keptvox)=p_FDR_mask;
            
            
            p_FWE_mask=p<(p_thr/(size(peakR_perm,1)+length(lags)))
            p_FWE_mask_3d=nan(voxn,1);
            p_FWE_mask_3d(keptvox)=p_FWE_mask;
            
            
            p_3d_nii=peakR;
            p_3d_nii.img=p_3d;
            save_nii(p_3d_nii,[expdir exp '/fmri/nii/wholeBrain/isfc_seed/isfc_seed_lagcorr_'  num2str(min(lags)) '-' num2str(max(lags))  '_' rname '_p.nii']);
            
            r_FWE=peakR;
            r_FWE.img(p_FWE_mask_3d==0)=NaN;
            t_FWE=peakT;
            t_FWE.img(p_FWE_mask_3d==0)=NaN;
            save_nii(t_FWE,[expdir exp '/fmri/nii/wholeBrain//aud/audenv_' subj '_'  num2str(min(lags)) '-' num2str(max(lags))  '_peakT_p' num2str(p_thr) 'FWE.nii']);
            save_nii(r_FWE,[expdir exp '/fmri/nii/wholeBrain//aud/audenv_' subj '_'  num2str(min(lags)) '-' num2str(max(lags))  '_peakR_p' num2str(p_thr) 'FWE.nii']);
            
            r_FDR=peakR;
            r_FDR.img(p_FDR_mask_3d==0)=NaN;
            t_FDR=peakT;
            t_FDR.img(p_FDR_mask_3d==0)=NaN;
            save_nii(t_FDR,[expdir exp '/fmri/nii/wholeBrain//aud/audenv_' subj '_'  num2str(min(lags)) '-' num2str(max(lags))  '_peakT_p' num2str(p_thr) 'FDR.nii']);
            save_nii(r_FDR,[expdir exp '/fmri/nii/wholeBrain//aud/audenv_' subj '_'  num2str(min(lags)) '-' num2str(max(lags))  '_peakR_p' num2str(p_thr) 'FDR.nii']);
            
        end
        
        clear peakR_perm peakT peakR
    end
end





