function roi2rois_bined_temporal_lagcorr_LL_leave1out(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
seed=rnames{ri};

crop_start=10;
lags_tested={-10:10, -40:40, -60:60};
binSize=250;

for ei=3:4;
    exp=experiments{ei};
    mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/perm/']);
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' seed '.mat' ],'gdata');
        gdata_seed=gdata;
        gdata_seed(:,:,subjects_excluded{ei})=NaN;
        gdata_seed=nanmean(gdata_seed,1);
        
        [~,tn,listenerN]=size(gdata_seed);
        
        r=nan([length(rnames)  tn length(lags) listenerN ]);
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
                gdata(:,:,subjects_excluded{ei})=NaN;
                gdata=nanmean(gdata,1);
                
                for si=1:listenerN;
                    
                    othersi=1:listenerN;
                    othersi=othersi(othersi~=si);
                    
                    for t=1:tn;
                        t_bin=t:(t+binSize-1);
                        
                        if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                            
                            y=nanmean(gdata(:,t_bin,othersi),3);
                            x=gdata_seed(:,t_bin,si);
                            
                            r(ri,t,:,si)=lagcorr(x',y',lags);
                        end
                    end
                end
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/' seed '_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames');
        clear r
    end
end


