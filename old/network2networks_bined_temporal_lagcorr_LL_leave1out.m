function network2networks_bined_temporal_lagcorr_LL_leave1out

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
seeds=networks;
crop_start=25;
crop_end=20;
lags_tested={-20:20, -40:40, -60:60};
binSize_tested=[30 20 100];

for ei=9:10;%[1 2 4 11 12];%[9 10];[3 5 7 8 9:12];
    exp=experiments{ei};
    mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/perm/']);
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        for bi=1;%1:length(binSize_tested);
            binSize=binSize_tested(bi);
            
            for sdi=1:length(networks);
                seed=networks{sdi};
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed '.mat' ],'gdata');
                gdata_seed=gdata;
                gdata_seed(:,:,subjects_excluded{ei})=NaN;
                gdata_seed=nanmean(gdata_seed,1);
                
                [~,tn,listenerN]=size(gdata_seed);
                keptT=(crop_start+1):(tn-crop_end);
                
                r=nan([length(networks) tn length(lags) listenerN ]);
                
                for ni=1:size(networks);
                    network=networks{ni};
                    
                    if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ]);
                        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ],'gdata');
                        gdata(:,:,subjects_excluded{ei})=NaN;
                        gdata=nanmean(gdata,1);
                        
                        for si=1:listenerN;
                            othersi=1:listenerN;
                            othersi=othersi(othersi~=si);
                            
                            for t=1:tn;
                                t_bin=t:(t+binSize-1);
                                
                                if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                                    
                                    y=nanmean(zscore(gdata(:,t_bin,othersi),0,2),3);
                                    x=gdata_seed(:,t_bin,si);
                                    
                                    r(ni,t,:,si)=lagcorr(x',y',lags);
                                    
                                end
                            end
                        end
                    end
                end
                save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/' seed '_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','networks','keptT');
                clear r
            end
        end
    end
end


