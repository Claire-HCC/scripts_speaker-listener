function roi2rois_temporal_lagcorr_LL_leave1out

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
seeds={'pANG_L','HG_L','vPCUN'};

crop_start=10;
lags_tested={-10:10, -40:40};

for ei=3:4;
    exp=experiments{ei};
%     rmdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/perm/'],'s');
%     mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/perm/']);
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        for sdi=1;%1:length(seed);
            seed=seeds{sdi};
            
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' seed '.mat' ],'gdata');
            gdata_seed=gdata;
            gdata_seed=nanmean(gdata_seed,1);
            [~,tn,listenerN]=size(gdata_seed);
            
           keptT_s=min(find(([1:tn]+min(lags))>0))+crop_start;
            keptT_e=min(tn,tn-max(lags));
            keptT=keptT_s:keptT_e;
            
            r=nan([length(rnames)  length(lags) listenerN ]);
            
            for ri=1:size(rnames);
                rname=rnames{ri};
                
                if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
                    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
                    gdata=nanmean(gdata,1);
                    
                    for si=1:listenerN;
                        othersi=1:listenerN;
                        othersi=othersi(othersi~=si);
                        
                        y=nanmean(nanmean(gdata(:,keptT,othersi),3),1);
                        x=gdata_seed(:,keptT,si);
                        
                        r(ri,:,si)=lagcorr(y',x',lags);
                        
                    end
                end
            end
            
            save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT');
            clear x y r
        end
    end
end

