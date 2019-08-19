clear all
close all
set_parameters

tic

% binary mask
load([expdir '/roi_mask/mat/merlin_isc_slm_pFDR01.mat'],'data');
mask=data;

for ei=1%:2;
    exp= experiments{ei};
    
    % 2d speaker data. vox X time
    load([expdir exp '/fmri/mat/wholeBrain/speaker01.mat']);
    speaker=zscore(data')';
    
    % 3d data. vox x time x listeners
    load([expdir experiments{ei} '/fmri/mat/wholeBrain/listenerAll.mat'],'data');
    listeners=data;
    clear data
    
    vox_selected=(mask==1 & sum(speaker')'~=0  & sum(sum(listeners==0,3),2)==0);
    
    % compute listener-speaker isfc for each listener
    for si = 1:size(listeners,3);
        % choose listener and take zscore of timeseries in each voxel.
        listener=zscore(listeners(:,:,si)')';
        % the output is a 3d mat. % listenr vox X speaker vox X listeners. 
        isfc_slm(:,:,si)= corr(listener(vox_selected,:)',speaker(vox_selected,:)');
    end
end
isfc=isfc_slm;
save([expdir exp '/fmri/mat/wholeBrain/isfc/isfc_slm_allSubj.mat' ],'isfc','vox_selected','-v7.3');

isfc=mean(isfc_slm,3);
save([expdir exp '/fmri/mat/wholeBrain/isfc/isfc_slm_mean.mat' ],'isfc','vox_selected','-v7.3');


