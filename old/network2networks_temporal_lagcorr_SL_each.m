function network2networks_temporal_lagcorr_SL_each

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
seeds=networks;
crop_start=10;
lags_tested={-10:10, -40:40, -60:60};

for ei=1:4;%1:4;
    exp=experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/SL_each/perm/']);
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        for sdi=1:length(networks);
            seed=networks{sdi};
            
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/speaker_' seed '.mat'],'data');
            data_seed=zscore(nanmean(data,1),0,2);
            
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed '.mat' ],'gdata');
            [~,tn,listenerN]=size(gdata);
            keptT=(crop_start+1):tn;
            
            r=nan([length(networks)  length(lags) listenerN ]);
            
            for ni=1:size(networks);
                network=networks{ni};
                
                if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ]);
                    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ],'gdata');
                    gdata(:,:,subjects_excluded{ei})=NaN;
                    gdata=zscore(nanmean(gdata,1),0,2);
                    
                    for si=1:listenerN;
                        
                        y=gdata(:,keptT,si);
                        x=data_seed(:,keptT);
                        
                        r(ni,:,si)=lagcorr(x',y',lags);
                    end
                end
            end
            save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/SL_each/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','networks','keptT');
            clear r
        end
    end
end


