clear all
close all
set_parameters

roles={'listenerZscoreMean','speaker01'}; % 'speaker'
itern=10000;

data=load_nii([expdir 'roi_mask/gray_matter_mask.nii']);
mask=data.img;

for rolei=1:length(roles);
    role=roles{rolei};
    
    for ei=1%:4;
        exp=experiments{ei};
        subjects=cellstr(ls([expdir exp '/fmri/mat/wholeBrain/' role '.mat']));
        % subjects= subjects(~cellfun(@isempty,(regexp(subjects,'[0-9]*'))));
        
        for si=1:length(subjects);
            
            % subj=subjects{si};
            [filepath,subj,ext] = fileparts(subjects{si});
            
            load([expdir experiments{ei} '/sound/' exp '_listener_cropped_audenv.mat' ]);
            aud=zscore(aud);
            
            load([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj '.mat']);
            epi=zscore(data')';
            
            keptvox=find(sum(data')'~=0 & mask(:)~=0);
            voln=min(size(epi,2),length(aud));
            
            %     peakR_perm=zeros(length(keptvox),length(lags),itern);
            peakR_perm=zeros(length(keptvox),1,itern);
            
            for iteri=1:itern;
                % phase_rand is done by column
                aud_perm= phase_rand(aud, 1);
                aud_perms=repmat(aud_perm(1:voln),1,length(keptvox))';
                %    lagcc_temp= lagcorr_claire(aud_perms',epi(keptvox,1:voln)',lags)';
                corr_temp= corr_col(aud_perms',epi(keptvox,1:voln)')';
                
                peakR_perm(:,:,iteri)=corr_temp;
                
                clear peakR
                
            end
        end
        
        save([expdir exp '/fmri/mat/wholeBrain/aud/audenv_' subj '_perm.mat'],'peakR_perm','keptvox','-v7.3');
        clear peakR_perm
    end
end



