
loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='restFc_isc30PercMasked_cluster6';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};

crop_start=25;
crop_end=20;
permN=3000;
for ei=1%:11;
    exp=exp_parameters.experiments{ei};
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/perm/']);
    
    for sdi=1;%2:length(networks);
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
                
                for perm=1:permN;
                    rng(perm)
                    
                    gdata_seed_perm=[];
                    gdata_perm=[];
                    gdata_seed_perm(1,:,:) =phase_rand2(squeeze(gdata_seed),0);
                    gdata_perm(1,:,:)=phase_rand2(squeeze(gdata),0);
                    
                    for si=1:listenerN;
                        
                        othersi=1:listenerN;
                        othersi=othersi(othersi~=si);
                        
                        y=nanmean(zscore(gdata_perm(:,keptT,othersi),0,2),3);
                        x=fliplr(zscore(gdata_seed_perm(:,keptT,si),0,2));
                        
                        r(tgi,:,si,perm)=circularlagcorr(x',y',lags);
                        
                    end
                end
            end
        end
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/perm/' seed '_permPhase'  ],'r','networks','keptT');
    end
end


