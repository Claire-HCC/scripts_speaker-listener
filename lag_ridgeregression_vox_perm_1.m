function []=lag_ridgeregression_vox_perm

clear all
close all

loc='cluster';
set_parameters;
load([expdir 'roi_mask/mat/MNI152_T1_3mm_brain_mask.mat']);
mask=data;
lags=[-10:10]; %scan

% params.iterations = 100;
params.iterations_PerPiece = 5;
params.iterations_piece = 1; % params.iterations/params.iterations_PerPiece;
params.permute      = 1;

for ei=1%:2;
    
    exp= experiments{ei};
    load([expdir exp '/fmri/mat/wholeBrain/speaker01.mat']);
    speaker=zscore(data')';
    
    % subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/listener*mat']));
    subjects=listeners;
    
    for itp=1:params.iterations_piece;
        
        B_sl=zeros(voxn,length(lags),params.iterations_PerPiece);
        B_ll=zeros(voxn,length(lags),params.iterations_PerPiece);
        
        for si=1:length(subjects);
            
            subj=subjects{si}
            
            load([expdir exp '/fmri/mat/wholeBrain/' subj ]);
            listener_self=zscore(data')';
            
            load([expdir exp '/fmri/mat/wholeBrain//leave1out/leave_' subj ]);
            listener_others=zscore(data')';
            
            for iti=1:params.iterations_PerPiece;
                params.iterations_PerPiece;
                
                disp('Running coupling filter and lag correlation analysis...');
                
                listener_self_perm = phase_rand(listener_self', params.permute)';
                
                vox_selected=100000;
               % vox_selected=(mask==1 & sum(speaker')'~=0  & sum(listener_self')'~=0);
                [B,F,p] = lag_ridgeregression3(listener_self_perm(vox_selected,:), speaker(vox_selected,:), [lags(1) lags(end)]);
                B_sl(vox_selected,:,iti)=B+ B_sl(vox_selected,:,iti);
                
                % vox_selected=(mask==1 &  sum(listener_others')'~=0  & sum(listener_self')'~=0);
                [B,F,p]= lag_ridgeregression3(listener_self_perm(vox_selected,:), listener_others(vox_selected,:), [lags(1) lags(end)]);
                B_ll(vox_selected,:,iti)=B+B_ll(vox_selected,:,iti);
            end
        end
        
        clear B F p
        B=B_sl/length(subjects);
        
        perfiles=ls([expdir experiments{ei} '/fmri/mat/wholeBrain/isc/perm/isc_sl_perm*']);
        iti=(regexp(perfiles,'[0-9]*','match'))'
        iti=cellfun(@str2num,iti);
        iti=max(iti);
        if isempty(iti); iti=0; end
        perm_idx=[num2str(iti+1) '-' num2str(iti+params.iterations_PerPiece) ];
        save([expdir experiments{ei} '/fmri/mat/wholeBrain/isc/perm/isc_sl_perm' perm_idx '.mat'],'B');
        disp('saving one piece...');
        clear B B_sl
        
        perfiles=ls([expdir experiments{ei} '/fmri/mat/wholeBrain/isc/perm/isc_ll_perm*']);
        iti=(regexp(perfiles,'[0-9]*','match'))'
        iti=cellfun(@str2num,iti);
        iti=max(iti);
             if isempty(iti); iti=0; end
        B=B_ll/length(subjects);
        save([expdir experiments{ei} '/fmri/mat/wholeBrain/isc/perm/isc_ll_perm' perm_idx '.mat'],'B');
        
        clear B B_ll 
    end
end

