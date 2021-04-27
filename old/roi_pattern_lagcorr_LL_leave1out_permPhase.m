function roi_pattern_lagcorr_LL_leave1out_permPhase(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
rname=rnames{ri};
crop_start=10;
lags_tested={-10:10, -40:40};
permN=1000;

for ei=3:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata')
        [~,tn,listenerN]=size(gdata);
        
        keptT=(crop_start+1):tn;
        
        r=nan([1  length(lags) listenerN permN]);
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            gdata(:,:,subjects_excluded{ei})=NaN;
            
            for perm=1:permN;
                rng(perm)
                gdata_perm=nan(size(gdata));
                for si=1:listenerN;
                    gdata_perm(:,keptT,si)=(phase_rand2(gdata(:,keptT,si)',1))';
                end
                
                for si=1:listenerN;
                    if ~(ei==1 & si==10); % that subject is excluded.
                        othersi=1:listenerN;
                        othersi=othersi(~ismember(othersi,si));
                        
                        y=nanmean(gdata_perm(:,keptT,othersi),3);
                        x=gdata_perm(:,keptT,si);
                        
                        [r(1,:,si,perm) ]=lagcorr_spatialTemporal(y,x,lags);
                    end
                end
            end
        end
        
        save(sprintf('%s/%s/fmri/pattern/lagcorr/%s/roi/%s/LL_leave1out/perm/lag%d-%d_permPhase_%s',expdir,exp,timeUnit,froidir,min(lags),max(lags),rname),'r','lags','rnames','keptT','-v7.3');
   
    end
end


