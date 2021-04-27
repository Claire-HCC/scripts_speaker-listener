

loc='mypc';
set_parameters;
timeUnit='tr' ;
crop_start=25;
crop_end=20;

froidir='restFc_isc30PercMasked_75Overlap_cluster6';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};

for ei=1%:8;
    exp=exp_parameters.experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/aud2networks//Ls/']);
    
    for sdi=1:length(networks);
        seed=networks{sdi};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed '.mat' ],'gdata','keptvox');
        
        gdata(:,:,exp_parameters.subjects_excluded{ei})=NaN;
        gdata=nanmean(gdata,1);
        
        [~,tn,listenerN]=size(gdata);
        keptT=(crop_start+1):(tn-crop_end);
        lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
        
        load([expdir exp '/sound/' exp '_listener_audhrf.mat' ],'aud');
        % I alreay cropped 3 TRs from the fMRI signal to account for the hrf
        % delay.
        aud=aud(4:(tn+3));
        
        if sdi==1;
            r=nan([length(networks) length(lags) listenerN ]);
        end
        for si=1:listenerN;
            y=gdata(:,keptT,si)';
            
            r(sdi,:,si)=circularlagcorr(aud(keptT),y,lags)';
        end
    end
    save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/aud2networks/Ls/aud.mat'   ],'r','lags','keptnetworks','keptT','networks');
end