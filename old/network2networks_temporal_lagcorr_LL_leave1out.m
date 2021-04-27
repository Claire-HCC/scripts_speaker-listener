function network2networks_temporal_lagcorr_LL_leave1out

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
seeds=networks;
crop_start=25;
crop_end=20;
lags_tested={-20:20, -15:15, -10:10, -40:40, -60:60};

for ei=[1 2 4 9:12];%1:12;%1:12;
    exp=experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks_peakLag0_mask/' froidir '/LL_leave1out/perm/']);
    
    for lagi=3;%:length(lags_tested);
        lags=lags_tested{lagi};
        
        for sdi=1:length(networks);
            seed=networks{sdi};
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed '_peakLag0_mask.mat' ],'gdata');
            gdata_seed=gdata;
            
          [~,tn,listenerN]=size(gdata_seed);
            keptT=(crop_start+1):(tn-crop_end);
            
            gdata_seed(:,:,subjects_excluded{ei})=NaN;
            gdata_seed=nanmean(gdata_seed,1);
            

            r=nan([length(networks)  length(lags) listenerN ]);
            
            for ni=1:size(networks);
                network=networks{ni};
                
                if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '_peakLag0_mask.mat' ]);
                    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '_peakLag0_mask.mat' ],'gdata');
            
                    gdata(:,:,subjects_excluded{ei})=NaN;
                    gdata=nanmean(gdata,1);
                    
                    for si=1:listenerN;
                        othersi=1:listenerN;
                        othersi=othersi(othersi~=si);
                        
                        y=nanmean(zscore(gdata(:,keptT,othersi),0,2),3);
                        x=zscore(gdata_seed(:,keptT,si),0,2);
                        
                        r(ni,:,si)=circularlagcorr(x',y',lags);
                        
                    end
                end
            end
            save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks_peakLag0_mask/' froidir '/LL_leave1out/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags))  ],'r','lags','networks','keptT');
            clear r
        end
    end
end


