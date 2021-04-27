function roi_pattern_lagcorr_LL_pairwise(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10, -30:30};


exp=experiments{ei};

for lagi=1:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rnames{1} '.mat' ],'gdata')
    [~,tn,listenerN]=size(gdata);
    
    keptT=(crop_start+1):tn;
    
    r=nan([length(rnames)  length(lags) listenerN listenerN]);
    
    for ri=1:size(rnames);
        rname=rnames{ri};
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
            
            for s1=1:listenerN;
                for s2=1:1:listenerN;
                    if ~(ei==1 & (s1==10 | s2==10)) & (s1~=s2) ; % that subject is excluded.
                        
                        y=zscore(gdata(:,keptT,s2),0,2);
                        x=zscore(gdata(:,keptT,s1),0,2);
                        
                        [r(ri,:,s1,s2) ]=lagcorr_spatialTemporal(y,x,lags);
                    end
                end
            end
        end
    end
    mkdir([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/LL_pairwise/perm/']);
    save([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/LL_pairwise/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT');
    clear x y r
end


