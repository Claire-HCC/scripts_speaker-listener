function network2networks_temporal_circularlagcorr_LL_leave1out_goodSubjs

loc='cluster';
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
for ei=[2 9:12];
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2vox/LL_leave1out/isc'   ],'peakLags','peaks','keptvox');
    maskPeakLag0=zeros(voxn,1);
    thr=sort(peaks(peakLags==0),'descend');
    thr=thr(round(length(keptvox)*0.3));
    maskPeakLag0(keptvox(peaks>thr & peakLags'==0))=1;
    clear peakLags peaks keptvox
    
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_goodSubjs/perm/']);
    
    load([expdir '/' exp '/bhv/comprehensionScore.mat'  ],'score');
    for sdi=1:length(networks);
        seed=networks{sdi};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed '.mat' ],'gdata','keptvox');
        % keep only voxel with peakLag 0
        disp(sum(ismember(keptvox,find(maskPeakLag0))))
        gdata_seed=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
        gdata_seed(:,:,subjects_excluded{ei})=NaN;
        gdata_seed(:,:,(score/max(score))<0.5)=NaN;
        gdata_seed=nanmean(gdata_seed,1);
        
        [~,tn,listenerN]=size(gdata_seed);
        keptT=(crop_start+1):(tn-crop_end);
        lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
        
        r=nan([length(networks)  length(lags) listenerN ]);
        
        for ni=1:size(networks);
            network=networks{ni};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ],'gdata','keptvox');
                % keep only voxels with peakLag 0
                gdata=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
                gdata(:,:,subjects_excluded{ei})=NaN;
                gdata(:,:,(score/max(score))<0.5)=NaN;
                gdata=nanmean(gdata,1);
                
                for si=1:listenerN;
                    othersi=1:listenerN;
                    othersi=othersi(othersi~=si);
                    
                    y=nanmean(zscore(gdata(:,keptT,othersi),0,2),3);
                    x=zscore(gdata_seed(:,keptT,si),0,2);
                    
                    r(ni,:,si)=circularlagcorr(x',y',lags);
                    
                end
            end
        end
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_goodSubjs/' seed   ],'r','networks','keptT');
    end
    
end


