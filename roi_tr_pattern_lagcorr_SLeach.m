function roi_tr_pattern_lagcorr_SLeach

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10, -30:30};

for ei=2:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rnames{1} '.mat' ],'gdata')
        [~,tn,listenerN]=size(gdata);
        
        keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
        keptT_e=min(tn,tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        r=nan([length(rnames)  length(lags) listenerN ]);
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rname '.mat'],'data');
                x=zscore(data(:,keptT),0,2);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
                roi_voxn=size(gdata,1);
                
                for si=1:listenerN;
                    y=gdata(:,keptT,si);
                    y=zscore(y,0,2);
                    
                    [r(ri,:,si) ]=lagcorr_spatialTemporal(y,x,lags);
                    
                    clear x y
                end
            end
        end
        % mkdir([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLeach/']);
        save([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLeach/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT');
        clear x y r
    end
end

