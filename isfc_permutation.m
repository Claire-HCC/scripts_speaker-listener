clear all
close all
set_parameters

tic

load([expdir '/roi_mask/mat/merlin_isc_slm_pFDR01.mat'],'data');
mask=data;

for ei=1%:2;
    exp= experiments{ei};
    load([expdir exp '/fmri/mat/wholeBrain/speaker01.mat']);
    speaker=zscore(data')';
    
    load([expdir experiments{ei} '/fmri/mat/wholeBrain/listenerAll.mat'],'data');
    listeners=data;
    for si = 1:size(listeners,3);
    listeners(:,:,si)=zscore(listeners(:,:,si)')';
    end

    subjn=size(listeners,3);
    
    group=[ones(1,subjn/2) 2*ones(1,subjn/2)];
    
    for iter=1:100;
        group=group(randperm(subjn));
        listener_g1=mean(listeners(:,:,group==1),3);
        listener_g2=mean(listeners(:,:,group==2),3);
        
        vox_selected=(mask==1 & sum(speaker')'~=0  & sum(listener_g1')'~=0 &  sum(listener_g2')'~=0);
      
        isfc_slm_g1= corr(listener_g1(vox_selected,:)',speaker(vox_selected,:)');
        isfc_slm_g2= corr(listener_g2(vox_selected,:)',speaker(vox_selected,:)');
        
        isfc_reliability(iter)=corr(isfc_slm_g1(:),isfc_slm_g2(:));
    end
end

save([expdir exp '/fmri/mat/wholeBrain/isfc/isfc_slm_reliability.mat' ],'isfc_reliability','-v7.3');


