clear all
close all
set_parameters

tic

% binary mask
load([expdir '/roi_mask/mat/merlin_isc_slm_pFDR01.mat'],'data');
mask=data;
lags= -10:10;
for ei=1%:2;
    exp= experiments{ei};
    
    % 2d speaker data. vox X time
    load([expdir exp '/fmri/mat/wholeBrain/speaker01.mat']);
    speaker=zscore(data')';
    
    % 3d data. vox x time x listeners
    load([expdir experiments{ei} '/fmri/mat/wholeBrain/listenerMean.mat'],'data');
    listeners=data;
    clear data
    
    vox_selected=find(mask==1 & sum(speaker')'~=0  & sum(sum(listeners==0,3),2)==0);
    
    % compute listener-speaker isfc for each listener
    for si = 1:size(listeners,3);
        % choose listener and take zscore of timeseries in each voxel.
        listener=zscore(listeners(:,:,si)')';
        
        for i=1:length(vox_selected);
            vi=vox_selected(i);
            % the output is a 3d mat. % listenr vox X speaker vox X listeners.
            listener_v=repmat(listener(vi,:)',1,length(vox_selected));
            
            % negative lag means listener precedes
            lagcc_temp= lagcorr_claire(speaker(vox_selected,:)',listener_v,lags);
            [peaksR peaksT]=max(lagcc_temp);
            peaksR_sorted=sort(peaksR,'descend');
            
            isfc_peakT(i,:,si)=mean(lags(peaksT(find(peaksR > peaksR_sorted(100)))));
            isfc_peakR(i,:,si)=mean(peaksR(find(peaksR > peaksR_sorted(100))));
        end
    end
end

save([expdir exp '/fmri/mat/wholeBrain/isfc/isfc_lagcorr_meanSubj_peakT.mat' ],'isfc_peakT','vox_selected','-v7.3');
save([expdir exp '/fmri/mat/wholeBrain/isfc/isfc_lagcorr_meanSubj_peakR.mat' ],'isfc_peakR','vox_selected','-v7.3');

toc
beep