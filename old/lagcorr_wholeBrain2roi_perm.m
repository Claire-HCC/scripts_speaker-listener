

% role1 roi to role2 wholebrain
clear all
close all
loc='mypc';
set_parameters

role1='listener01';
role2='listener01_others';
relation='LL_s01_leave1out';
lags=[-10:10]; %scan
itern=1000;

rnames={'HG_L','PMC_L','vPCUN','vmPFC','aANG_L','STC_L','pIFG_L'};

mask=  load_nii([expdir 'roi_mask/rspm8_gray0.5_merlin_groupAvtive_mask.nii']);
mask=mask.img(:);

for ei=1%:2;%1:2;
    exp=experiments{ei};
    
    load([expdir experiments{ei} '/fmri/mat/wholeBrain/' role1 '.mat']);
    data1=data;
    
    for ri=3%4:length(rnames);
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
        
        keptvox=find(mask==1 & sum(sum(data2==0,3),2)==0);
        
        peakR_perm=zeros(length(keptvox),1,itern);
        % peakR_perm=zeros(length(keptvox),length(lags),itern);
        
        for iteri=1:itern;
            % phase_rand is done by column
            role1_roi_perm= phase_rand(role1_roi', 1)';
            role1_roi_perms=repmat(role1_roi_perm,length(keptvox),1);
            
            %     lagcc_temp= lagcorr_claire(role1_roi_perms',listener(keptvox,:)',lags)';
            corr_temp= corr_col(role1_roi_perms',data2(keptvox,:)')';
            
            peakR_perm(:,:,iteri)=corr_temp;
            clear peakR role1_roi_perm role1_roi_perms corr_temp
        end
        
        save([expdir experiments{ei} '/fmri/mat/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_' rname '_perm.mat'],'peakR_perm','keptvox','-v7.3');
        
    end
    % end
end


beep