function roi2rois_bined_temporal_lagcorr_LL_leave1out_permPhase(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
seed='vPCUN';

crop_start=10;
lags_tested={-10:10, -40:40, -60:60};
binSize=30;
permN=1000;
rname=rnames{ri};

for ei=1:4;
    exp=experiments{ei};
    mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/perm/']);
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' seed '.mat' ],'gdata');
        gdata_seed=gdata;
        gdata_seed(:,:,subjects_excluded{ei})=NaN;
        gdata_seed=nanmean(gdata_seed,1);
        
        [~,tn,listenerN]=size(gdata_seed);
        
        r=nan([1  tn  length(lags) listenerN  permN ]);
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            gdata(:,:,subjects_excluded{ei})=NaN;
            gdata=nanmean(gdata,1);
            
            for perm=1:permN;
                
                gdata_seed_perm=nan(size(gdata_seed));
                gdata_perm=nan(size(gdata));
                for si=1:listenerN;
                    gdata_seed_perm(1,(crop_start+1):tn,si)=(phase_rand2(gdata_seed(:,(crop_start+1):tn,si)',1))';
                    gdata_perm(1,(crop_start+1):tn,si)=(phase_rand2(gdata(:,(crop_start+1):tn,si)',1))';
                end
                
                for si=1:listenerN;
                    othersi=1:listenerN;
                    othersi=othersi(othersi~=si);
                    
                    for t=1:tn;
                        t_bin=t:(t+binSize-1);
                        
                        if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                            
                            y=nanmean(gdata_perm(:,t_bin,othersi),3);
                            x=gdata_seed_perm(:,t_bin,si);
                            
                            r(1,t,:,si,perm)=lagcorr(y',x',lags);
                        end
                    end
                end
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/perm/' seed '_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_' rname],'r','lags','rnames');
        clear r
    end
end


