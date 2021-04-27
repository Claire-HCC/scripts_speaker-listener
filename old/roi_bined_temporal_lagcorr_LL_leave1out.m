function roi2roi_temporal_lagcorr_LL_leave1out

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10, -40:40, -60:60};
binSize=20:

for ei=1;%1:4;
    exp=experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2roi/' froidir '/LL_leave1out/perm/']);
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata');
        
        [~,tn,listenerN]=size(gdata);
        
        r=nan([length(rnames) tn length(lags) listenerN ]);
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
                gdata(:,:,subjects_excluded{ei})=NaN;
                gdata=nanmean(gdata,1);
                
                for si=1:listenerN;
                    for t=1:tn;
                        t_bin=t:(t+binSize-1);
                        
                        if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                            
                            othersi=1:listenerN;
                            othersi=othersi(othersi~=si);
                            
                            y=nanmean(gdata(:,t_bin,othersi),3);
                            x=gdata(:,t_bin,si);
                            
                            r(ri,t,:,si)=lagcorr(y',x',lags);
                        end
                    end
                end
            end
        end
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2roi/' froidir '/LL_leave1out_bined/bin' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT');
        clear r
    end
end


