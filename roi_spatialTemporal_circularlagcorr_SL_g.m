function roi_spatialTemporal_circularlagcorr_SL_g
% canonical correlation is always high. probably because of too many voxels
% in the ROI mask, compared to the limited number of time points we have.
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
crop_start=25;
crop_end=20;

for ei=1%:2;%[1 2];%[1 2 4 9:12]
    exp=experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/pattern/circularlagcorr/' timeUnit '/roi/' froidir '/SL_g/']);
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rnames{1} '.mat' ],'data','keptvox');
    [~,tn]=size(data);
    keptT=(crop_start+1):(tn-crop_end);
    lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
    
    r=nan([length(rnames)  length(lags)  ]);
    
    for ri=1:length(rnames);
        
        rname=rnames{ri};
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
            
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rname '.mat' ],'data','keptvox');
            data_seed=data;
            
            [~,tn]=size(data);
            keptT=(crop_start+1):(tn-crop_end);
            lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
            
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata','keptvox');
            
            y=nanmean(zscore(gdata(:,keptT,:),0,2),3);
            x=zscore(data_seed(:,keptT),0,2);
            
            r(ri,:)= circularlagcorr_spatialTemporal(x,y,lags);
            
        end
    end
end
save([expdir '/' exp '/fmri/pattern/circularlagcorr/' timeUnit '/roi/' froidir '/SL_g/r.mat'],'r','lags','rnames','keptT');
clear r
end



