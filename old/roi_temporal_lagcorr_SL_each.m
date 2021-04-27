function roi_temporal_lagcorr_SL_each

 loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10, -40:40};

for ei=1:4;
    exp=experiments{ei};
%     rmdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_each/'],'s');
%     mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_each/perm/']);
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rnames{1} '.mat' ],'gdata')
        [~,tn,listenerN]=size(gdata);
        
        keptT=(crop_start+1):tn;
        
        r=nan([length(rnames)  length(lags) listenerN]);
        p=nan([length(rnames)  length(lags)]);
        t=p;
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rname '.mat'],'data');
                gdata(:,:,subjects_excluded{ei})=NaN;
                data=nanmean(data);
                gdata=nanmean(gdata);
                
                for si=1:listenerN;
                    y=gdata(:,keptT,si);
                    x=data(:,keptT);
                    
                    [r(ri,:,si)]=lagcorr(x',y',lags);
                    
                end
                [~,p(ri,:),~,stats]=ttest(squeeze(r(ri,:,:))');
                t(ri,:)=stats.tstat;
            end
        end
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT','p','t');
        clear r
    end
end
