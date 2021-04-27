function network2networks_temporal_lagcorr_LL_pairwise

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
seeds=networks;
crop_start=20;
crop_end=20;
lags_tested={-10:10, -40:40, -60:60};

for ei=9;
    exp=experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_pairwise/perm/']);
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        for sdi=1:length(networks);
            seed=networks{sdi};
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed '.mat' ],'gdata');
            gdata_seed=gdata;
            gdata_seed(:,:,subjects_excluded{ei})=NaN;
            gdata_seed=nanmean(gdata_seed,1);
            
            [~,tn,listenerN]=size(gdata_seed);
            
            keptT=(crop_start+1):(tn-crop_end);
            
            r=nan([length(networks)  length(lags) listenerN listenerN]);
            
            for ni=1:size(networks);
                network=networks{ni};
                
                if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ]);
                    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ],'gdata');
                    gdata(:,:,subjects_excluded{ei})=NaN;
                    gdata=nanmean(gdata,1);
                    
                    for si1=1:listenerN;
                        for si2=1:listenerN;
                            if si1~=si2;
                                y=gdata(:,keptT,si2);
                                x=gdata_seed(:,keptT,si1);
                                
                                r(ni,:,si1,si2)=lagcorr(x',y',lags);
                                
                            end
                        end
                    end
                end
            end
            save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_pairwise/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','networks','keptT');
            clear r
        end
    end
end


