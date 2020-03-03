function roi2rois_tr_temporal_lagcorr_LL_leave1out_permPhase(ri)

%loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
seed='HG_L';
crop_start=10;
lags_tested={-10:10, -40:40};
permN=1000;
rname=rnames{ri};

for ei=3;
    exp=experiments{ei};
    
    for lagi=2;%:2%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' seed '.mat' ],'gdata');
        gdata_seed=gdata;
        gdata_seed=nanmean(gdata_seed,1);
        
        [~,tn,listenerN]=size(gdata_seed);
        
        keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
        keptT_e=min(tn,tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        r=nan([1  length(lags) listenerN permN]);
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            gdata=nanmean(gdata,1);
            
            for perm=1:permN;
                rng(perm)
                
                seed_perm=[];
                for si=1:listenerN;
                    if si==1;
                        [temp rp]=(phase_rand2(gdata_seed(:,:,si)',1));
                        seed_perm(:,:,si)=temp';
                    else
                        seed_perm(:,:,si)=(phase_rand2(gdata_seed(:,:,si)',1,rp))';
                    end
                end
                for si=1:listenerN;
                    othersi=1:listenerN;
                    othersi=othersi(othersi~=si);
                    y=nanmean(gdata(:,keptT,othersi),3);
                    
                    x=seed_perm(:,keptT,si);
                    r(1,:,si,perm)=lagcorr(y',x',lags);
                end
            end
        end
        
        %     mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/perm/']);
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/perm/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permPhase1_' rname],'r','lags','keptT');
    end
end

