function roi_temporal_circularlagcorr_SL_g
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'],'data');
    
loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10, -40:40};

for ei=1:2;
    exp=experiments{ei};
    rmdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_gg/perm/'],'s');
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi/' froidir '/SL_g/perm/']);
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_gg/perm/']);
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rnames{1} '.mat' ],'data')
        tn=size(data,2);
        
        keptT=(crop_start+1):tn;
        
        r=nan([length(rnames)  length(lags) ]);
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            
                gdata(:,:,subjects_excluded{ei})=NaN;
                
                y=nanmean(nanmean(gdata(:,keptT,:),3));
                x=nanmean(data(:,keptT));
                
                [r(ri,:) ]=circularlagcorr(y,x,lags);
            end
        end
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT');
    end
end

