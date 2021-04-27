function roi2wholeBrain_temporal_lagcorr_LL_leave1out

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

crop_start=10;
lags_tested={-10:10, -30:30};
seed='vPCUN';

for ei=4;%[1:2]
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' seed '.mat' ],'gdata','keptvox');
        gdata_seed=gdata;
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll_zscore.mat' ],'gdata','keptvox');
        [voxn,tn,listenerN]=size(gdata);
             keptT=(crop_start+1):tn;
        
        r=nan([voxn length(lags) listenerN]);
        
        for si=1:listenerN;
            othersi=1:listenerN;
            othersi=othersi(othersi~=si);
            
            y=nanmean(gdata(:,keptT,othersi),3);
            x=repmat(mean(gdata_seed(:,keptT,si),1),length(keptvox),1);
            
            r(:,:,si)=(lagcorr_claire(y',x',lags))';
        end
        mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2wholeBrain/LL_leave1out/perm']);
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2wholeBrain/LL_leave1out/'  seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','keptvox','keptT','-v7.3');
        % clear r
    end
end


