function roi_temporal_lagcorr_SL_g_permPhase(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
rname=rnames{ri};

crop_start=10;
lags_tested={-10:10,  -40:40};
permN=1000;

for ei=1:2;%2:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'],'data');
            gdata(:,:,subjects_excluded{ei})=NaN;
            
            gdata=nanmean(gdata,1);
            data=nanmean(data,1);
            
            [~,tn,~]=size(gdata);
            keptT=(crop_start+1):tn;
            
            y=nanmean(gdata(:,keptT,:),3);
            
            r=nan([1  length(lags) permN]);
            for perm=1:permN;
                x=phase_rand2(data(:,keptT)',1)';
                [r(1,:,perm) ,~]=lagcorr(y',x',lags);
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_' rname],'r','lags','keptT','-v7.3');
        clear r
    end
end

