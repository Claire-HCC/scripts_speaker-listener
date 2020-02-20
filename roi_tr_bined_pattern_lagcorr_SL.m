clear all;

tic
loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
binSize_tested=[40]; % tr;
lags_tested={-10:10 , -30:30};

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

% cropt start because there is clearly a speech-start effect in the
% listeners' data

for ei=1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rnames{1} '.mat'],'data');
    tn=size(data,2);
    
    for binSizei=1:length(binSize_tested);
        binSize=binSize_tested(binSizei);
        
        for lagi=2;%1:length(lags_tested);
            lags=lags_tested{lagi};
            
            r=nan([length(rnames) tn  length(lags) ]);
            
            
            for ri=1:length(rnames);
                rname=rnames{ri};
                
                if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ]);
                    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
                    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
                    
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
            mkdir([expdir '/' exp '/fmri/pattern_lagcorr_bined/' timeUnit '/roi/' froidir '/SLg/perm/']);
            save([expdir '/' exp '/fmri/pattern_lagcorr_bined/' timeUnit '/roi/' froidir '/SLg/SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','rnames','binSize','lags');
            clear r
        end
    end
end
toc
beep
