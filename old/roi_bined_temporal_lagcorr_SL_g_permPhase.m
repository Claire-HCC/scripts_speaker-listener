function roi_bined_temporal_lagcorr_SL_g_permPhase(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
binSize_tested=[30]; % tr;
lags_tested={-10:10 , -40:40};

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
rname=rnames{ri};

crop_start=10;
% cropt start because there is clearly a speech-start effect in the
% listeners' data
permN=1000;
for ei=1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rnames{1} '.mat'],'data');
    tn=size(data,2);
    permN=1000;
    
    for binSizei=1%:length(binSize_tested);
        binSize=binSize_tested(binSizei);
        
        for lagi=1%;1:length(lags_tested);
            lags=lags_tested{lagi};
            
            r=nan([1 tn  length(lags) permN ]);
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'],'data');
                gdata(:,:,subjects_excluded{ei})=NaN;
                gdata=nanmean(gdata);
                data=nanmean(data);
                g=nanmean(gdata,3);
                
                for perm=1:permN;
                    data_perm=nan(size(data));
                    data_perm(:,(crop_start+1):tn)=(phase_rand2(data(:,(crop_start+1):tn)',1))';
                    
                    for t=1:tn;
                        t_bin=t:(t+binSize-1);
                        
                        if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                            
                            y=g(:,t_bin);
                            x=data_perm(:,t_bin);
                            [r(1,t,:,perm) ]=lagcorr(y',x',lags);
                        end
                    end
                end
            end
            save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_g_bined/perm/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_' rname ],'r','rnames','binSize','lags');
            clear r
        end
    end
end


