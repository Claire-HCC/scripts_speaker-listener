function roi_temporal_lagcorr_SL_g_permSL(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10, -40:40};

exp=experiments{ei};
%  mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi/' froidir '/SL_g/perm/']);

for lagi=1%:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata');
    [~,tn,listenerN]=size(gdata);
    r=nan([length(rnames)  length(lags) listenerN ]);
    keptT=(crop_start+1):tn;
    
    for ri=1:size(rnames);
        rname=rnames{ri};
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'],'data');
            gdata(:,:,subjects_excluded{ei})=NaN;
            roi_voxn=size(gdata,1);
            
            for perm=1:listenerN;
                data_perm=gdata(:,:,perm);
                gdata_perm=gdata;
                gdata_perm(:,:,perm)=data;
                
                g=nanmean(gdata_perm,3);
                y=g(:,keptT);
                
                x=data_perm(:,keptT);
                
                [r(ri,:,perm) ]=lagcorr(y',x',lags);
                
                clear x y
                
            end
        end
    end
    save( sprintf('%s/%s/fmri/temporal/lagcorr/%s/roi/%s/SL_g/perm/lag%d-%d_permSL' , expdir,exp, timeUnit,froidir,min(lags),max(lags)),'r','lags','rnames','keptT');
end



