function roi2rois_tr_temporal_lagcorr_SS

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
seed='vPCUN';
crop_start=10;
lags_tested={-10:10, -40:40};

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' seed '.mat' ],'data');
        data_seed=data;
        data_seed=nanmean(data_seed,1);
        [~,tn]=size(data_seed);
        
        keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
        keptT_e=min(tn,tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        r=nan([length(rnames)  length(lags) ]);
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat' ],'data');
                gdata=nanmean(gdata,1);
                
                othersi=1:listenerN;
                othersi=othersi(othersi~=si);
                
                y=nanmean(nanmean(gdata(:,keptT,othersi),3),1);
                x=gdata_seed(:,keptT,si);
                
                r(ri,:)=lagcorr(y',x',lags);
            end
            
        end
        mkdir([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi2rois/' froidir '/SS/']);
        save([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi2rois/' froidir '/SS/' seed 'lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','rnames','keptT');
        clear x y r
    end
end

