function roi_tr_pattern_lagcorr_LL_leave1out_permPhase(perm)
rng(perm)
loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10, -30:30};

for ei=[3 4 1 2]
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rnames{1} '.mat' ],'gdata')
        [~,tn,listenerN]=size(gdata);
        
        keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
        keptT_e=min(tn,tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        r=nan([length(rnames)  length(lags) listenerN]);
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
                
                for si=1:listenerN;
                    if ~(ei==1 & si==10); % that subject is excluded.
                        othersi=1:listenerN;
                        othersi=othersi(~ismember(othersi,si));
                        
                        roi_voxn=size(gdata,1);
                        
                        g=nanmean(gdata(:,:,othersi),3);
                        y=g;
                        y=zscore(y(:,keptT),0,2);
                        
                        x=phase_rand2(gdata(:,:,si)',1);
                        x=x';
                        x=zscore(x(:,keptT),0,2);
                        
                        [r(ri,:,si) ]=lagcorr_spatialTemporal(y,x,lags);
                    end
                    clear x y
                end
            end
        end
        
        save(sprintf('%s/%s/fmri/pattern_lagcorr/%s/roi/%s/LL_leave1out/perm/lag%d-%d_permPhase%04d',expdir,exp,timeUnit,froidir,min(lags),max(lags),perm),'r','lags','rnames','keptT');
        clear x y r
    end
end


