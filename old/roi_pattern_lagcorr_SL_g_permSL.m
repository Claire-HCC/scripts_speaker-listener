function roi_temporal_lagcorr_SL_g_permSL(perm)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10};

for ei=3:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata');
        [~,tn,listenerN]=size(gdata);
        
        if perm<=listenerN;
            keptT=(crop_start+1):tn;
            
            r=nan([length(rnames)  length(lags) ]);
            
            for ri=1:size(rnames);
                rname=rnames{ri};
                
                if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
                    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
                    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'],'data');
                    
                    roi_voxn=size(gdata,1);
                    
                    data_perm=gdata(:,:,perm);
                    gdata_perm=gdata;
                    gdata_perm(:,:,perm)=data;
                    
                    g=nanmean(gdata_perm,3);
                    y=g;
                    y=zscore(y(:,keptT),0,2);
                    
                    x=zscore(data_perm(:,keptT),0,2);
                    
                    [r(ri,:) ]=lagcorr_spatialTemporal(y,x,lags);
                    
                    clear x y
                    
                end
            end
            
            save( sprintf('%s/%s/fmri/temporal/lagcorr/%s/roi/%s/SL_g/perm/lag%d-%d_permSL%02d' , expdir,exp, timeUnit,froidir,min(lags),max(lags),perm),'r','lags','rnames','keptT');
            clear x y r
        end
    end
end

