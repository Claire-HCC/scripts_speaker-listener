function roi_tr_temporal_lagcorr_SLg_permPhase(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
rname=rnames{ri};

crop_start=10;
lags_tested={-10:10,  -30:30};
permN=1000;

for ei=1;%1:4%2:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/permPhase/speaker_' rname '.mat'],'data');
            data_perm=data;
            
            [~,tn,~]=size(gdata);
            keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
            keptT_e=min(tn,tn-max(lags));
            keptT=keptT_s:keptT_e;
            
            gdata=nanmean(gdata,1);
            g=nanmean(gdata,3);
            y=g;
            y=y(:,keptT);
            
            r=nan([1  length(lags) permN]);
            for perm=1:permN;
                data=data_perm(:,:,perm);
                % average all voxels within the roi
                data=nanmean(data,1);
                
                x=data(:,keptT);
                
                [r(1,:,perm) ,~]=lagcorr(y',x',lags);
                
                clear x
            end
        end
        
        save([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/perm/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_' rname],'r','lags','keptT');
        clear x y r
    end
end

