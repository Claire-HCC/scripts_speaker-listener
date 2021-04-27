function roi2rois_temporal_circularlagcorr_LL_gg

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
seeds={'HG_L','vPCUN','pANG_L'};

crop_start=10;
lags_tested={-10:10, -40:40, -60:60};

for ei=1:4;
    exp=experiments{ei};
    %  rmdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_gg/perm/'],'s');
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_gg/perm/']);
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        for sdi=1:length(seeds);
            seed=seeds{sdi};
            
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' seed '.mat' ],'gdata');
            gdata_seed=gdata;
            gdata_seed(:,:,subjects_excluded{ei})=NaN;
            
            [~,tn,~]=size(gdata_seed);
            keptT=(crop_start+1):tn;
            
            gdata_seed=nanmean(gdata_seed(:,keptT,:),1);
            
            r=nan([length(rnames)  length(lags) ]);
            
            for ri=1:size(rnames);
                rname=rnames{ri};
                
                if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]) & strcmpi(rname,seed)==0;
                    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
                    gdata(:,:,subjects_excluded{ei})=NaN;
                    gdata=nanmean(gdata(:,keptT,:),1);
                    
                    y=nanmean(gdata,3);
                    x=nanmean(gdata_seed,3);
                    
                    r(ri,:)=circularlagcorr(y',x',lags);
                    
                end
            end
            save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_gg/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT');
            clear r
        end
    end
end


