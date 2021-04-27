function roi2rois_temporal_lagcorr_SL_each_permPhase(sdi)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
seeds=rnames;
permN=1000;
crop_start=10;
lags_tested={-10:10, -40:40, -60:60};

for ei=[3 4 ]
    exp=experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/SL_each/perm/']);
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        %  for sdi=1:length(seeds);
        seed=seeds{sdi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' seed '.mat' ],'gdata');
        gdata_seed=gdata;
        gdata_seed(:,:,subjects_excluded{ei})=NaN;
        gdata_seed=nanmean(gdata_seed,1);
        
        [~,tn,listenerN]=size(gdata_seed);
        keptT=(crop_start+1):tn;
        
        r=nan([length(rnames)  length(lags) listenerN ]);
        
          rng(perm)
                    
                    gdata_seed_perm=[];
                    gdata_perm=[];
                    for si=1:listenerN;
                        gdata_seed_perm(1,:,si)=(phase_rand2(gdata_seed(:,:,si)',1))';
                        gdata_perm(1,:,si)=(phase_rand2(gdata(:,:,si)',1))';
                    end
                    
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
                gdata(:,:,subjects_excluded{ei})=NaN;
                gdata=nanmean(gdata,1);
                
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rname '.mat' ],'data');
                data=nanmean(data,1);
                
                for si=1:listenerN;
                    
                    y=gdata(:,keptT,si);
                    x=data(:,keptT);
                    
                    r(ri,:,si)=lagcorr(y',x',lags);
                    
                end
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/SL_each/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT');
        clear r
        % end
    end
end

