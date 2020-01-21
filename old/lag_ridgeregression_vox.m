clear all
close all
set_parameters
load([expdir 'roi_mask\mat\MNI152_T1_3mm_brain_mask.mat']);
mask=data;

lags=[-10:10]; %scan


for ei=1%:2;
    
    exp= experiments{ei};
    load([expdir exp '/fmri/mat/wholeBrain/speaker01.mat']);
    speaker=zscore(data')';
    
    subjects=listeners;
    
    %     B_sl=zeros(voxn,length(lags),length(subjects));
    %     F_sl=zeros(voxn,length(subjects));
    %     p_sl=zeros(voxn,length(subjects));
    
    B_slm=zeros(voxn,length(lags));
    F_slm=zeros(voxn,1);
    p_slm=zeros(voxn,1);
    
    %     B_ll=zeros(voxn,length(lags),length(subjects));
    %     F_ll=zeros(voxn,length(subjects));
    %     p_ll=zeros(voxn,length(subjects));
    
    %     for si=1:length(subjects);
    %         subj=subjects{si};
    %
    %         load([expdir exp '/fmri/mat/wholeBrain/' subj ]);
    %         listener_self=zscore(data')';
    %
    %         load([expdir exp '/fmri/mat/wholeBrain//leave1out/leave_' subj ]);
    %         listener_others=zscore(data')';
    %
    %         disp('Running coupling filter and lag correlation analysis...');
    %
    %         vox_selected=(mask==1 & sum(speaker')'~=0  & sum(listener_self')'~=0);
    %
    %         [B,F,p] = lag_ridgeregression3(listener_self(vox_selected,:), speaker(vox_selected,:), [lags(1) lags(end)]);
    %         B_sl(vox_selected,:,si)=B;
    %         F_sl(vox_selected,si)=F;
    %         p_sl(vox_selected,si)=p;
    %
    %         vox_selected=(mask==1 &  sum(listener_others')'~=0  & sum(listener_self')'~=0);
    %         [B,F,p]= lag_ridgeregression3(listener_self(vox_selected,:), listener_others(vox_selected,:), [lags(1) lags(end)]);
    %         B_ll(vox_selected,:,si)=B;
    %         F_ll(vox_selected,si)=F;
    %         p_ll(vox_selected,si)=p;
    %
    %     end
    
    load([expdir exp '/fmri/mat/wholeBrain/listenerZscoreMean.mat' ]);
    listenerMean=zscore(data')';
    
    vox_selected=(mask==1 & sum(speaker')'~=0  & sum(listenerMean')'~=0);
    
    [B,F,p] = lag_ridgeregression3(listenerMean(vox_selected,:), speaker(vox_selected,:), [lags(1) lags(end)]);
    B_slm(vox_selected,:)=B;
    F_slm(vox_selected)=F;
    p_slm(vox_selected)=p;
    
    
    clear B F p
    
    %     B=B_sl; F=F_sl; p=p_sl;
    %     save([expdir experiments{ei} '/fmri/mat/wholeBrain/isc/isc_sl.mat'],'B','F','p');
    %     B=B_ll; F=F_ll; p=p_ll;
    %     save([expdir experiments{ei} '/fmri/mat/wholeBrain/isc/isc_ll.mat'],'B','F','p');
    
    B=B_slm; F=F_slm; p=p_slm;
    save([expdir experiments{ei} '/fmri/mat/wholeBrain/isc/isc_slm.mat'],'B','F','p','vox_selected');
end

nii=mat2nii(F);
save_nii(nii,[expdir experiments{ei} '/fmri/nii/wholeBrain/isc/isc_slm_F.nii']);

vox_selected=(F~=0);
p_temp=sort(p(vox_selected))
FDR = mafdr(p_temp);
fdri=max(find(FDR<=0.01));
p_FDR=p_temp(fdri);
F_thr=finv(1-p_FDR,21,585-22)

data=F>F_thr;
save([expdir '/roi_mask/mat/merlin_isc_slm_pFDR01.mat'],'data');
nii=mat2nii(data);
save_nii(nii,[expdir  '/roi_mask/nifti/merlin_isc_slm_pFDR01.nii']);

