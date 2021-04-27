function roi_pattern_circularlagcorr_SL_g

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10, -40:40};

for ei=1:4;
    exp=experiments{ei};
    mkdir([expdir '/' exp '/fmri/pattern/circularlagcorr/' timeUnit '/roi/' froidir '/SL_g/perm/']);
    
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
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'],'data');
                gdata(:,:,subjects_excluded{ei})=NaN;
                
                y=nanmean(gdata(:,keptT,:),3);
                x=data(:,keptT);
                
                [r(ri,:) ]=circularlagcorr_spatialTemporal(y,x,lags);
            end
        end
        
        save([expdir '/' exp '/fmri/pattern/circularlagcorr/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT');
    end
end

