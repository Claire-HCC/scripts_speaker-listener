function roi_tr_bined_pattern_lagcorr_SL_permPhase(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
binSize_tested=[40]; % tr;
lags_tested={-10:10 };

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
rname=rnames{ri}
% cropt start because there is clearly a speech-start effect in the
% listeners' data
permN=1000;


for ei=1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rnames{1} '.mat'],'data');
    tn=size(data,2);
    
    for binSizei=1:length(binSize_tested);
        binSize=binSize_tested(binSizei);
        
        for lagi=1:length(lags_tested);
            lags=lags_tested{lagi};
            
            r=nan([1 tn  length(lags) permN]);
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/permPhase/speaker_' rname '.mat'],'data');
                data=zscore(data,0,2);
                
                roi_voxn=size(gdata,1);
                g=nanmean(gdata(:,:,:),3);
                
                for perm=1:permN;
                    
                    for t=1:tn;
                        t_bin=t:(t+binSize-1);
                        
                        if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                            % substract the global mean pattern
                            y=g(:,t_bin);
                            x=data(:,t_bin,perm);
                            
                            [r(1,t,:,perm) ]=lagcorr_spatialTemporal(y,x,lags);
                            
                            clear x y
                        end
                    end
                end
                
             
                save([expdir '/' exp '/fmri/pattern_lagcorr_bined/' timeUnit '/roi/' froidir '/SLg/perm/SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_' rname],'r','rnames','binSize','lags','-v7.3');
                clear r
            end
        end
    end
end
