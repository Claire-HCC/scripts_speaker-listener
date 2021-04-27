function network2networks_temporal_circularlagcanoncorr_SLg

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
seeds=networks;
crop_start=25;
crop_end=20;

fsize=[35 18];
figure('unit','centimeter','position',[0 0 fsize]);
for ei=[1 2];%[4 11 9 10];
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2vox/LL_leave1out/isc'   ],'peakLags','peaks','keptvox');
    maskPeakLag0=zeros(voxn,1);
    thr=sort(peaks(peakLags==0),'descend');
    thr=thr(round(length(keptvox)*0.3));
    maskPeakLag0(keptvox(peaks>thr & peakLags'==0))=1;
    clear peakLags peaks keptvox
    
    
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcanoncorr/' timeUnit '/network2networks/' froidir '/SLg/']);
    
    for sdi=11;%:length(networks);
        seed=networks{sdi};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/speaker_' seed '.mat' ],'data','keptvox');
        % keep only voxel with peakLag 0
        disp(sum(ismember(keptvox,find(maskPeakLag0))))
        gdata_seed=data(ismember(keptvox,find(maskPeakLag0)),:,:);
        gdata_seed(:,:,subjects_excluded{ei})=NaN;
        
        [~,tn,listenerN]=size(gdata_seed);
        keptT=(crop_start+1):(tn-crop_end);
        lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
        
        r=nan([ length(networks) length(lags)  ]);
        
        for ni=1:size(networks);
            network=networks{ni};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ],'gdata','keptvox');
                % keep only voxels with peakLag 0
                gdata=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
                gdata(:,:,subjects_excluded{ei})=NaN;
                
                y=nanmean(zscore(gdata(:,keptT,:),0,2),3);
                x=zscore(gdata_seed(:,keptT),0,2);
                
                r(ni,:)=circularlagcorr_spatialTemporal_canonical(x,y,lags);
                
            end
        end
        save([expdir '/' exp '/fmri/temporal/circularlagcanoncorr/' timeUnit '/network2networks/' froidir '/SLg/' seed   ],'r','networks','keptT');
    end
    
end




