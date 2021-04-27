function roi_temporal_lagcorr_LL_leave1out_permPhase(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
crop_start=10;
lags_tested={-10:10, -40:40};
permN=1000;
rname=rnames{ri};

for ei=[3 4 ]
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata');
        
        [~,tn,listenerN]=size(gdata);
        keptT=(crop_start+1):tn;
        
        r=nan([1  length(lags) listenerN permN]);
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            gdata(:,:,subjects_excluded{ei})=NaN;
            gdata=nanmean(gdata,1);
            
            for perm=1:permN;
                rng(perm)
                
                gdata_perm=nan(size(gdata));
                for si=1:listenerN;
                    gdata_perm(1,keptT,si)=(phase_rand2(gdata(:,keptT,si)',1))';
                end
                
                for si=1:listenerN;
                    othersi=1:listenerN;
                    othersi=othersi(othersi~=si);
                    
                    y=nanmean(gdata_perm(:,:,othersi),3);
                    x=gdata_perm(:,:,si);
                    
                    r(1,:,si,perm)=lagcorr(y',x',lags);
                end
            end
        end
        mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/LL_leave1out/perm']);
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/LL_leave1out/perm/lag' num2str(min(lags)) '-' num2str(max(lags))  '_permPhase_' rname],'r','lags','keptT');
        clear r
    end
end



