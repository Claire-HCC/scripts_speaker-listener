function roi2rois_temporal_lagcorr_SL_g_permPhase(sdi)

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
crop_start=10;
lags_tested={-10:10, -40:40, -60:60};
seeds={'vPCUN','HG_L','pANG_L'};
seed=seeds{sdi};
permN=2;

for ei=[3 4 ]
    exp=experiments{ei};
    %    mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/SL_g/perm/']);
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' seed '.mat' ],'data');
        [~,tn]=size(data);
        keptT=(crop_start+1):tn;
        
        x=zscore(nanmean(data(:,keptT)),0,2);
        r=nan([length(rnames)  length(lags)  permN]);
        
        for perm=1:permN;
            rng(perm)
            
            x_perm=phase_rand2(x',1)';
            
            for ri=1%:size(rnames);
                rname=rnames{ri};
                
                if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
                    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
                    gdata(:,:,subjects_excluded{ei})=NaN;
                    gdata=nanmean(gdata,1);
                    
                    y=nanmean(zscore(gdata(:,keptT,:),0,2),3);
                    
                    r(ri,:,perm)=lagcorr(y',x_perm',lags);
                    
                end
            end
        end
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/SL_g/perm/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags))   '_permPhase'],'r','lags','rnames','keptT');
        clear r
    end
end

