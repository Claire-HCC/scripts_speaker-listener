function roi2rois_temporal_circularlagcorr_SL_g(sdi)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
seeds=rnames;
seed=seeds{sdi};
crop_start=25;
crop_end=20;

for ei=[1 2];%[1 2 4 9:12]
    exp=experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/SL_g/']);
    
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2vox/LL_leave1out/isc'   ],'peakLags','peaks','keptvox');
    maskPeakLag0=zeros(voxn,1);
    thr=sort(peaks(peakLags==0),'descend');
    thr=thr(round(length(keptvox)*0.3));
    maskPeakLag0(keptvox(peaks>thr & peakLags'==0))=1;
    clear peakLags peaks keptvox
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' seed '.mat' ],'data','keptvox');
    
    [~,tn]=size(data);
    keptT=(crop_start+1):(tn-crop_end);
    lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
    r=nan([length(rnames)  length(lags)  ]);
    
    % keep only voxel with peakLag 0
    if (sum(ismember(keptvox,find(maskPeakLag0)))) > 10;
        data_seed=data(ismember(keptvox,find(maskPeakLag0)),:,:);
        data_seed=nanmean(data_seed,1);
        
        for ri=1:size(rnames);
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata','keptvox');
                if (sum(ismember(keptvox,find(maskPeakLag0)))) > 10;
                    gdata=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
                    gdata(:,:,subjects_excluded{ei})=NaN;
                    gdata=nanmean(gdata,1);
                    
                    y=nanmean(zscore(gdata(:,keptT,:),0,2),3);
                    x=zscore(data_seed(:,keptT),0,2);
                    
                    r(ri,:)=circularlagcorr(x',y',lags);
                    
                end
            end
        end
    end
    save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/SL_g/' seed ],'r','lags','rnames','keptT');
    clear r
end



