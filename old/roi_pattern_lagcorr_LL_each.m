function roi_pattern_lagcorr_LL_each(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10, -40:40};

exp=experiments{ei};
mkdir(sprintf('%s/%s/fmri/pattern/lagcorr/%s/roi/%s/LL_each/', expdir,exp, timeUnit,froidir));

for lagi=1%:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata');
    [~,tn,listenerN]=size(gdata);
    r=nan([length(rnames)  length(lags) listenerN listenerN ]);
    keptT=(crop_start+1):tn;
    
    for ri=1:size(rnames);
        rname=rnames{ri};
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            gdata(:,:,subjects_excluded{ei})=NaN;
            roi_voxn=size(gdata,1);
            
            for si1=1:listenerN;
                for si2=1:listenerN;
                    if si1~=si2;
                        y=gdata(:,keptT,si2);
                        x=gdata(:,keptT,si1);
                        [r(ri,:,si1,si2) ]=lagcorr_spatialTemporal(y,x,lags);
                        clear x y
                    end
                end
            end
        end
    end
    save( sprintf('%s/%s/fmri/pattern/lagcorr/%s/roi/%s/LL_each/lag%d-%d' , expdir,exp, timeUnit,froidir,min(lags),max(lags)),'r','lags','rnames','keptT');
end



