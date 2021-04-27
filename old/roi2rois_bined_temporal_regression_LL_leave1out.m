function roi2rois_bined_temporal_regression_LL_leave1out(ei)

 loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
seeds={'vPCUN','HG_L'};

lags_tested={-10:10, -40:40};
binSize=30;

exp=experiments{ei};
mkdir([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/perm/']);

for lagi=1%:length(lags_tested);
    lags=lags_tested{lagi};
    
    for sdi=1%:length(seeds);
        seed=seeds{sdi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' seed '.mat' ],'gdata');
        gdata_seed=gdata;
        gdata_seed(:,:,subjects_excluded{ei})=NaN;
        gdata_seed=nanmean(gdata_seed,1);
        
        [~,tn,listenerN]=size(gdata_seed);
        roi_voxn=1;
        
        b=nan([length(rnames) tn length(lags)+1 listenerN]);
        r2=nan([length(rnames) tn listenerN]);
        F=nan([length(rnames) tn listenerN]);
        p=nan([length(rnames) tn listenerN]);
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
                gdata(:,:,subjects_excluded{ei})=NaN;
                gdata=nanmean(gdata,1);
                
                
                for si=1:listenerN;
                    
                    othersi=1:listenerN;
                    othersi=othersi(othersi~=si);
                    
                    for t=1:tn;
                        t_bin=t:(t+binSize-1);
                        
                        if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                            
                            y=(nanmean(gdata(:,t_bin,othersi),3));
                            y=y(:);
                            for li=1:length(lags);
                                X(:,:,li)=gdata_seed(:,t_bin+lags(li),si);
                            end
                            
                            X=reshape(X,roi_voxn*length(t_bin),length(lags));
                            
                            % centralized X
                            X=X-mean(X);
                            
                            [b(ri,t,:,si),~,r,~,stats]=regress(y,[ones(size(X,1),1) X]);
                            
                            r2(ri,t,si)=stats(1);
                            F(ri,t,si)=stats(2);
                            p(ri,t,si)=stats(3);
                            
                            clear X
                        end
                    end
                end
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/' seed '_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','rnames','-v7.3');
        clear r
    end
end


