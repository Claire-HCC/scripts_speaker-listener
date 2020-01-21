close all
clear
loc='mypc';
set_parameters

roles={'listener','speaker01'}; % 'speaker'

mask_nii=load_nii([expdir '/roi_mask/gray_matter_mask.nii']);
mask=mask_nii.img;

for rolei=1%:length(roles);
    role=roles{rolei};
    
    for ei=1;
        exp=experiments{ei};
        subjects=cellstr(ls([expdir exp '/fmri/mat/wholeBrain/' role '*.mat']));
        subjects= subjects(~cellfun(@isempty,(regexp(subjects,'[0-9][0-9]'))));
        
        load([expdir exp '\sound\listenerEvents.mat']);
        X=[ones(size(listenerEvents_designMat,1),1) listenerEvents_designMat];
        
        tic
        for si=2:length(subjects);
            
            [filepath,subj,ext] = fileparts(subjects{si});
            modeldir=[expdir experiments{ei} '/fmri/nii/wholeBrain/L1/' subj ];
            mkdir(modeldir);
            
            load([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj '.mat']);
            epi=zscore(data,0,2);
            
            keptvox=find(mask(:)==1 & (sum(data~=0,2)==size(data,2)));
            data=data(keptvox,:);
            
            bs=[];
            betas=zeros(voxn,size(X,2));
            
            for kvi=1:length(keptvox);
                vi=keptvox(kvi);
                y=data(kvi,:)';
                if sum(y)~=0;
                    b=glmfit(X,y,'normal','constant','off')';
                    betas(vi,:)=b;
                end
                
            end
            save([modeldir '/betas.mat'],'betas','X')
            
            for bi=2:size(X,2);
                nii=mat2nii(betas(:,bi),loc);
                fname=sprintf('%s/beta_Event%02d.nii',modeldir,bi-1)
                save_nii(nii,fname);
            end
            clear data keptvox betas
        end
        
        toc
        beep
    end
end