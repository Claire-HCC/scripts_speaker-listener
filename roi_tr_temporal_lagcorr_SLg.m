function roi_tr_temporal_lagcorr_SLg

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-10:10, -30:30};

for ei=1%:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rnames{1} '.mat' ],'data')
        tn=size(data,2);
                
        keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
        keptT_e=min(tn,tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        r=nan([length(rnames)  length(lags) ]);
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rname '.mat'],'data');
                
                % average all voxels within the roi
                gdata=nanmean(gdata,1);
                data=nanmean(data,1);
                
                roi_voxn=size(gdata,1);
                
                g=nanmean(gdata,3);
                y=g;
                y=y(:,keptT);
                
                x=data(:,keptT);
                
                [r(ri,:) ,~]=lagcorr(y',x',lags);
                
                clear x y
                
            end
        end
        mkdir([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/']);
        save([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT');
        clear x y r
    end
end

