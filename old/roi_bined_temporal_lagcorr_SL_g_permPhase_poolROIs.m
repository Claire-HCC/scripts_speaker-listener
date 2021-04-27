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
            
            rtemp=nan([length(rnames) tn  length(lags) permN ]);
            for ri=1:length(rnames);
                rname=rnames{ri};
                load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_g_bined/perm/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_' rname ],'r','rnames','binSize','lags');
                rtemp(ri,:,:,:)=r;
            end
            r=rtemp;
            
            save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_g_bined/perm/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase'],'r','rnames','binSize','lags','-v7.3');
            
        end
    end
end

