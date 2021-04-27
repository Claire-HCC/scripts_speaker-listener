function network2networks_temporal_circularlagcorr_LL_selfself_timeReversed

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='restFc_isc30PercMasked_75Overlap_cluster6';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};
crop_start=25;
crop_end=20;

for ei=10;
    exp=exp_parameters.experiments{ei};
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_selfself/perm/']);
    
    for sdi=1:length(networks);
        seed=networks{sdi};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed '.mat' ],'gdata','keptvox');
        gdata_seed=gdata;
        gdata_seed(:,:,exp_parameters.subjects_excluded{ei})=NaN;
        gdata_seed=nanmean(gdata_seed,1);
        
        [~,tn,listenerN]=size(gdata_seed);
        keptT=(crop_start+1):(tn-crop_end);
        lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
        
        r=nan([length(networks)  length(lags) listenerN ]);
        
        for tgi=1:length(networks);
            network=networks{tgi};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ],'gdata','keptvox');
                gdata(:,:,exp_parameters.subjects_excluded{ei})=NaN;
                gdata=nanmean(gdata,1);
                
                for si=1:listenerN;
                    
                    y=zscore(gdata(:,keptT,si),0,2);
                    x=fliplr(zscore(gdata_seed(:,keptT,si),0,2));
                    
                    r(tgi,:,si)=circularlagcorr(x',y',lags);
                    
                end
            end
        end
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_selfself/' seed '_timeReversed'  ],'r','networks','keptT');
    end
    
end


