function roi_temporal_lagcorr_LL_leave1out

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10, -40:40, -60:60};

for ei=1:4;
    exp=experiments{ei};
%    rmdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/LL_leave1out/'],'s');
%         mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/LL_leave1out/perm/']);
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rnames{1} '.mat' ],'gdata');
        
        [~,tn,listenerN]=size(gdata);
        keptT=(crop_start+1):tn;
        
        r=nan([length(rnames)  length(lags) listenerN ]);
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
                gdata(:,:,subjects_excluded{ei})=NaN;
                gdata=nanmean(gdata,1);
                
                for si=1:listenerN;
                    othersi=1:listenerN;
                    othersi=othersi(othersi~=si);
                    
                    y=nanmean(gdata(:,keptT,othersi),3);
                    x=gdata(:,keptT,si);
                    
                    r(ri,:,si)=lagcorr(x',y',lags);
                    
                end
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/LL_leave1out/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT');
        clear r
    end
end


