function roi2rois_temporal_lagcorr_LL_leave1out(sdi)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
seeds=rnames;
seed=seeds{sdi};
crop_start=25;
crop_end=20;
lags_tested={-20:20, -40:40, -60:60};

for ei=[12];
    exp=experiments{ei};
    %  rmdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/'],'s');
    mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/perm/']);
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' seed '.mat' ],'gdata');
        gdata_seed=gdata;
        gdata_seed(:,:,subjects_excluded{ei})=NaN;
        gdata_seed=nanmean(gdata_seed,1);
        
        [~,tn,listenerN]=size(gdata_seed);
        keptT=(crop_start+1):(tn-crop_end);
        
        r=nan([length(rnames)  length(lags) listenerN ]);
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
                gdata(:,:,subjects_excluded{ei})=NaN;
                gdata=nanmean(gdata,1);
                
                for si=1:listenerN;
                    othersi=1:listenerN;
                    othersi=othersi(othersi~=si);
                    
                    y=nanmean(zscore(gdata(:,keptT,othersi),0,2),3);
                    x=zscore(gdata_seed(:,keptT,si),0,2);
                    
                    r(ri,:,si)=lagcorr(x',y',lags);
                    
                end
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT');
        clear r
        % end
    end
end

