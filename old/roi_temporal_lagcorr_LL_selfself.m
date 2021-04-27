function roi_temporal_lagcorr_LL_selfself(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10, -40:40};

exp=experiments{ei};
mkdir(sprintf('%s/%s/fmri/temporal/lagcorr/%s/roi/%s/LL_each/', expdir,exp, timeUnit,froidir));

for lagi=1%:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata');
    [~,tn,listenerN]=size(gdata);
    r=nan([length(rnames)  length(lags) listenerN  ]);
    keptT=(crop_start+1):tn;
    
    for ri=1:size(rnames);
        rname=rnames{ri};
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            gdata(:,:,subjects_excluded{ei})=NaN;
            gdata=nanmean(gdata,1);
            roi_voxn=1;
            
            for si=1:listenerN;
                x=gdata(:,keptT,si);
                [r(ri,:,si) ]=lagcorr(x',x',lags);
                clear x
            end
        end
    end
    mkdir(sprintf('%s/%s/fmri/temporal/lagcorr/%s/roi/%s/LL_selfself/' , expdir,exp, timeUnit,froidir))
    save( sprintf('%s/%s/fmri/temporal/lagcorr/%s/roi/%s/LL_selfself/lag%d-%d' , expdir,exp, timeUnit,froidir,min(lags),max(lags)),'r','lags','rnames','keptT');
end



