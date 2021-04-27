

loc='cluster';
set_parameters;
timeUnit='tr' ;
crop_start=25;
crop_end=20;

for ei=1:8;
    exp=exp_parameters.experiments{ei};
    
   
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/aud2vox//Ls/']);
    
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll_isc30PercMasked.mat',expdir,exp,timeUnit);
    load(f,'gdata','keptvox');
    gdata(:,:,exp_parameters.subjects_excluded{ei})=NaN;
    
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
    
    load([expdir exp '/sound/' exp '_listener_audhrf.mat' ],'aud');
    % I alreay cropped 3 TRs from the fMRI signal to account for the hrf
    % delay.
    aud=aud(4:(tn+3));
    
    r=nan(length(keptvox),length(lags),listenerN);
    for si=1:listenerN;
        y=nanmean(zscore(gdata(:,keptT,si),0,2),3)';
        x=repmat(aud(keptT),1,size(y,2));
        r(:,:,si)=circularlagcorr_byCol(x,y,lags)';
    end
    save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/aud2vox/Ls/aud.mat'   ],'r','lags','keptvox','keptT','-v7.3');
end
