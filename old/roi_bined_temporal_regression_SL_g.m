function roi_bined_temporal_regression_SL_g(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

lags_tested={-10:10, -40:40};
binSize=30;

exp=experiments{ei};
mkdir([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_g_bined/perm/']);

for lagi=1%:length(lags_tested);
    lags=lags_tested{lagi};
    
     load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata');
    [~,tn,listenerN]=size(gdata);
    roi_voxn=1;
    
    b=nan([length(rnames) tn length(lags)+1 );
    r2=nan([length(rnames) tn ]);
    F=nan([length(rnames) tn ]);
    p=nan([length(rnames) tn ]);
    
    for ri=1:size(rnames);
        rname=rnames{ri};
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat' ],'data');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            gdata(:,:,subjects_excluded{ei})=NaN;
            gdata=nanmean(gdata,1);
            data=nanmean(data);
            
            for t=1:tn;
                t_bin=t:(t+binSize-1);
                
                if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                    
                    y=(nanmean(gdata(:,t_bin,:),3));
                    y=y(:);
                    for li=1:length(lags);
                        X(:,:,li)=data(:,t_bin+lags(li));
                    end
                    
                    X=reshape(X,roi_voxn*length(t_bin),length(lags));
                    
                    % centralized X
                    X=X-mean(X);
                    
                    [b(ri,t,:),~,r,~,stats]=regress(y,[ones(size(X,1),1) X]);
                    
                    r2(ri,t)=stats(1);
                    F(ri,t)=stats(2);
                    p(ri,t)=stats(3);
                    
                    clear X
                end
            end
        end
    end
end

save([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_g_bined/' seed '_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','rnames','-v7.3');
clear r
end



