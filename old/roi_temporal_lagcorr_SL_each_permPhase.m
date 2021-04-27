function roi_temporal_lagcorr_SL_each_permPhase(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
rname=rnames{ri};
crop_start=10;
lags_tested={-10:10, -40:40};
permN=1000;

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata')
        [~,tn,listenerN]=size(gdata);
        
        r=nan([1  length(lags) listenerN permN]);
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'],'data');
            % gdata(:,:,subjects_excluded{ei})=NaN;
            data=nanmean(data,1);
            gdata=nanmean(gdata,1);
            
            roi_voxn=1;
            tn=size(data,2);
            
            keptT=(crop_start+1):tn;
            
            r=nan([length(rnames)  length(lags) listenerN]);
            for perm=1:permN;
                data_perm=nan(size(data));
                data_perm(:,keptT)=(phase_rand2(data(:,keptT)',1))';
                x=data_perm(:,keptT);
                
                for si=1:listenerN;
                    y=gdata(:,keptT,si);
                    
                    [r(1,:,si,perm)]=lagcorr(y',x',lags);
                    
                end
                
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_each/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_' rname],'r','lags','rnames','keptT');
        clear b F p r2 r2_byTime
    end
end


