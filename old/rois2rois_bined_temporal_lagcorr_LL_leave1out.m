function rois2rois_bined_temporal_lagcorr_LL_leave1out(ei)

 loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10, -40:40, -60:60};
binSize=30;

exp=experiments{ei};
mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/rois2rois/' froidir '/LL_leave1out_bined/perm/']);

for lagi=1%:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata');
    [~,tn,listenerN]=size(gdata);
    
    r=nan([length(rnames) length(rnames) tn length(lags) listenerN ]);
    
    for ri1=2%:length(rnames);
        rname1=rnames{ri1};
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname1 '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname1 '.mat' ],'gdata');
            gdata1=gdata;
            gdata1(:,:,subjects_excluded{ei})=NaN;
            gdata1=nanmean(gdata1,1);
            
            for ri2=11%1:size(rnames);
                rname2=rnames{ri2};
                
                if  exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname2 '.mat' ]);
                    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname2 '.mat' ],'gdata');
                    gdata2=gdata;
                    gdata2(:,:,subjects_excluded{ei})=NaN;
                    gdata2=nanmean(gdata2,1);
                    
                    for si=1:2;%1:listenerN;
                        
                        othersi=1:listenerN;
                        othersi=othersi(othersi~=si);
                        
                        for t=1:tn;
                            t_bin=t:(t+binSize-1);
                            
                            if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                                
                                y=nanmean(gdata2(:,t_bin,othersi),3);
                                x=gdata1(:,t_bin,si);
                                
                                r(ri1,ri2,t,:,si)=lagcorr(y',x',lags);
                            end
                        end
                    end
                end
            end
        end
    end
    save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/rois2rois/' froidir '/LL_leave1out_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','-v7.3');
    clear r
end



