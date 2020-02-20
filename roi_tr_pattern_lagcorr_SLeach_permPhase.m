function roi_tr_pattern_lagcorr_SLeach_permPhase(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10, -30:30};
permN=1000;
for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rnames{1} '.mat' ],'gdata')
        [~,tn,listenerN]=size(gdata);
        
        keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
        keptT_e=min(tn,tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        r=nan([1 length(lags) listenerN pernN]);
        
        rname=rnames{ri};
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/permPhase/speaker_' rname '.mat'],'data');
            
            for perm=1:permN;
                x=zscore(data(:,keptT,perm),0,2);
                
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
                roi_voxn=size(gdata,1);
                
                for si=1:listenerN;
                    y=gdata(:,keptT,si);
                    y=zscore(y,0,2);
                    
                    [r(1,:,si,perm) ]=lagcorr_spatialTemporal(y,x,lags);
                end
            end
        end
        rmdir([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLeach/perm/'],'s')
        mkdir([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLeach/perm/']);
        save(  sprintf('%s/%s/fmri/pattern_lagcorr/%s/roi/%s/SLeach/perm/SL_lag%d-%d_permPhase_%s',expdir,exp,timeUnit,froidir,min(lags),max(lags),rname),'r','lags','keptT');
        clear x y r
    end
end

