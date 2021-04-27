function roi_bined_pattern_lagcorr_SL(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
binSize_tested=[30]; % tr;
lags_tested={-10:10 , -40:40};

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

% cropt start because there is clearly a speech-start effect in the
% listeners' data


exp=experiments{ei};
mkdir([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_g_bined/perm/']);

load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rnames{1} '.mat'],'data');
tn=size(data,2);

for binSizei=1%:length(binSize_tested);
    binSize=binSize_tested(binSizei);
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        r=nan([length(rnames) tn  length(lags) ]);
        
        for ri=1:length(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'],'data');
                gdata(:,:,subjects_excluded{ei})=NaN;
                roi_voxn=size(gdata,1);
                g=nanmean(gdata(:,:,:),3);
                
                for t=1:tn;
                    t_bin=t:(t+binSize-1);
                    
                    if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                        % substract the global mean pattern
                        y=g(:,t_bin);
                        x=data(:,t_bin);
                        
                        [r(ri,t,:) ]=lagcorr_spatialTemporal(y,x,lags);
                        
                        clear x y
                    end
                end
            end
        end
        
        save([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_g_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','rnames','binSize','lags');
        clear r
    end
end


