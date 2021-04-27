function roi2rois_temporal_circularlagcorr_LL_leave1out_permPhase(ri)

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
seeds={'vPCUN','HG_L''pANG_L'};

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        for sdi=1:length(seeds);
            seed=seeds{sdi};
            
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' seed '.mat' ],'gdata');
            gdata_seed=gdata;
            
            [~,tn,listenerN]=size(gdata_seed);
            keptT=(crop_start+1):tn;
            
            gdata_seed(:,:,subjects_excluded{ei})=NaN;
            gdata_seed=nanmean(gdata_seed(:,keptT,:),1);
            
            r=nan([1  length(lags) listenerN permN]);
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
                
                gdata(:,:,subjects_excluded{ei})=NaN;
                gdata=nanmean(gdata(:,keptT,:),1);
                
                for perm=1%:permN;
                    rng(perm)
                    
                    gdata_seed_perm=[];
                    gdata_perm=[];
                    for si=1:listenerN;
                        gdata_seed_perm(1,:,si)=(phase_rand2(gdata_seed(:,:,si)',1))';
                        gdata_perm(1,:,si)=(phase_rand2(gdata(:,:,si)',1))';
                    end
                    
                    for si=1:listenerN;
                        othersi=1:listenerN;
                        othersi=othersi(othersi~=si);
                        
                        y=nanmean(gdata_perm(:,:,othersi),3);
                        x=gdata_seed_perm(:,:,si);
                        
                        r(1,:,si,perm)=circularlagcorr(y',x',lags);
                    end
                end
            end
            
            save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/perm/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permPhase_' rname],'r','lags','keptT');
            clear r
        end
    end
end


