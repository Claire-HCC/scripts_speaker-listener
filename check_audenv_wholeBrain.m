clear all
close all
loc='cluster';
set_parameters

roles={'listenerZscoreMean','speaker01'}; % 'speaker'
lags=-10:10;

data=load_nii([expdir 'roi_mask/gray_matter_mask.nii']);
mask=data.img;

for rolei=1:length(roles);
    role=roles{rolei};
    
    for ei=1;
        exp=experiments{ei};
        subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/' role '.mat']));
        % subjects= subjects(~cellfun(@isempty,(regexp(subjects,'[0-9]*'))));
        
        for si=1%:length(subjects);
            
            [filepath,subj,ext] = fileparts(subjects{si});
            
            %  load([expdir experiments{ei} '/sound/' exp '_' strrep(subj,'.mat','') '_audenv.mat' ]);
            load([expdir experiments{ei} '/sound/' exp '_listener_cropped_audenv.mat' ]);
            aud=zscore(aud);
            
            load([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj '.mat']);
            epi=zscore(data')';
            
            voln=min(length(aud),size(epi,2));
            
            keptvox=find(sum(data')'~=0 & mask(:)~=0);
            
            aud_mat=repmat(aud(1:voln),1,length(keptvox))';
            
            lagcc_temp= lagcorr_claire(aud_mat',epi(keptvox,1:voln)',lags)';
            [peakR peakT]=max(lagcc_temp');
            peakT=lags(peakT);
            
            save([expdir exp '/fmri/mat/wholeBrain/aud/audenv_lagcorr_' subj '_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'peakT','peakR','keptvox');
            
            peakR_nii=nan(voxn,1);
            peakR_nii(keptvox,:)=peakR;
            nii=mat2nii(peakR_nii);
            save_nii(nii,[expdir experiments{ei} '/fmri/nii/wholeBrain/aud/audenv_lagcorr_' subj '_' num2str(min(lags)) '-' num2str(max(lags)) '_peakR.nii']);
            
            peakT_nii=nan(voxn,1);
            peakT_nii(keptvox,:)=peakT;
            nii=mat2nii(peakT_nii);
            save_nii(nii,[expdir experiments{ei} '/fmri/nii/wholeBrain/aud/audenv_lagcorr_' subj '_' num2str(min(lags)) '-' num2str(max(lags)) '_peakT.nii']);
               
            clear peakR peakT
        end
    end
end


