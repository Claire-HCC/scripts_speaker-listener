function roi_pattern_circularlagcorr_SL_g_permPhase(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
rname=rnames{ri};
crop_start=10;
lags_tested={-10:10,  -40:40};
permN=1000;

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=2%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rnames{1} '.mat' ],'data')
        tn=size(data,2);
        
              keptT=(crop_start+1):tn;
        
        r=nan([1  length(lags) permN ]);
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]) & exist([expdir '/' exp '/fmri/pattern/circularlagcorr/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_' rname '.mat'])==0;
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'],'data');
            
            y=nanmean(gdata(:,keptT,:),3);
            
            for perm=1:permN;
                rng(perm);
                
                data_perm=(phase_rand(data',0))';
                x=data_perm(:,keptT);
                [r(1,:,perm) ]=circularlagcorr_spatialTemporal(y,x,lags);
            end
        end
        outf= sprintf('%s/%s/fmri/pattern/circularlagcorr/%s/roi/%s/SL_g/perm/lag%d-%d_permPhase_%s',expdir,exp,timeUnit,froidir,min(lags),max(lags),rname);
        save(outf,'r','lags','keptT','-v7.3');
        clear x y r
    end
end

