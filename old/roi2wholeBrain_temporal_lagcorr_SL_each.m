function roi2wholeBrain_temporal_lagcorr_SL_each

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

crop_start=10;
lags_tested={-10:10};
seed='vPCUN';

for ei=3;%
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' seed '.mat' ],'data')
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll.mat' ],'gdata','keptvox');
        [voxn,tn,listenerN]=size(gdata);
        r=nan(voxn,length(lags),listenerN);
        keptT=(crop_start+1):tn;
        
        for si=1:listenerN;
            x=repmat(mean(data(:,keptT),1),length(keptvox),1);
            y=gdata(:,keptT,si);
            
            r(:,:,si)=[lagcorr_claire(y',x',lags)]';
        end
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2wholeBrain/SL_each/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','keptvox','keptT','-v7.3');
    end
    
    clear x y r
end


