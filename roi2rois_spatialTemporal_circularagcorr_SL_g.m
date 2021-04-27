function roi2rois_spatialTemporal_circularcanonlagcorr_SL_g(sdi)
% canonical correlation is always high. probably because of too many voxels
% in the ROI mask, compared to the limited number of time points we have. 
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
seeds=rnames;
seed=seeds{sdi};
crop_start=25;
crop_end=20;

for ei=[1 2];%[1 2 4 9:12]
    exp=experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/pattern/circularlagcorr/' timeUnit '/roi2rois/' froidir '/SL_g/']);
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' seed '.mat' ],'data','keptvox');
    data_seed=data;
    
    [~,tn]=size(data);
    keptT=(crop_start+1):(tn-crop_end);
    lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
    r=nan([length(rnames)  length(lags)  ]);
    
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata','keptvox');
                
                y=nanmean(zscore(gdata(:,keptT,:),0,2),3);
                x=zscore(data_seed(:,keptT),0,2);
                
                r(ri,:)= circularlagcorr_spatialTemporal_canonical(x,y,lags);
                
            end
        end
    end
    save([expdir '/' exp '/fmri/pattern/circularcanonlagcorr/' timeUnit '/roi2rois/' froidir '/SL_g/' seed ],'r','lags','rnames','keptT');
    clear r
end



