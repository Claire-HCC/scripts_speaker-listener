function roi2rois_temporal_circularlagcorr_LL_gg_permPhase(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
seeds={'vPCUN','HG_L','pANG_L'};
rname=rnames{ri};
crop_start=10;
lags_tested={-10:10, -40:40};
permN=1000;

for ei=1:4;
    exp=experiments{ei};
    %  rmdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_gg/perm/'],'s');
    %  mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_gg/perm/']);
    
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
            
            r=nan([1  length(lags) permN]);
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
                gdata(:,:,subjects_excluded{ei})=NaN;
                gdata=nanmean(gdata(:,keptT,:),1);
                
                y=nanmean(gdata,3);
                x=nanmean(gdata_seed,3);
                
                for perm=1:permN;
                    x=phase_rand(x',1)';
                    r(1,:,perm)=circularlagcorr(y',x',lags);
                end
            end
            
            save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_gg/perm/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_' rname],'r','lags','rnames','keptT');
            clear r
        end
    end
end


