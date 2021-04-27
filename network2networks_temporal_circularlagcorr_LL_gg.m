function network2networks_temporal_circularlagcorr_LL_gg

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='fc_cluster';
% load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
% networks=unique(table2array(roi_table(:,2)));
networks={'network1','network2','network3','network4','network5','network6'};
seeds=networks;
crop_start=25;
crop_end=20;

iters=1000;

for ei=[7 8]
    exp=exp_parameters.experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_gg/perm/']);
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' networks{1} '.mat' ],'gdata');
    [~,tn,listenerN]=size(gdata);
     keptT=(crop_start+1):(tn-crop_end);
    lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
    
    subjs_g1=[];
    subjs_g2=[];
    for iter=1:iters;
        rng(iter)
        subjs_shuffled=randperm(listenerN);
        subjs_shuffled(ismember(subjs_shuffled,exp_parameters.subjects_excluded{ei}))=[];
        subjs_g1(:,iter)=subjs_shuffled(1:round(length(subjs_shuffled)/2));
        subjs_g2(:,iter)=subjs_shuffled((1+round(length(subjs_shuffled)/2)):end);
    end
    
    for sdi=1:length(networks);
        seed=networks{sdi};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed '.mat' ],'gdata');
        gdata_seed=gdata;
        gdata_seed=nanmean(gdata_seed,1);
        gdata_seed(:,:,exp_parameters.subjects_excluded{ei})=NaN;
        [~,tn,listenerN]=size(gdata_seed);
        
        keptT=(crop_start+1):(tn-crop_end);
        
        r=nan([length(networks)  length(lags) iters]);
        
        
        for ni=1:length(networks);
            network=networks{ni};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ],'gdata');
                gdata(:,:,exp_parameters.subjects_excluded{ei})=NaN;
                gdata=nanmean(gdata,1);
                
                for iter=1:iters;
                    x=nanmean(zscore(gdata_seed(:,keptT,subjs_g1(:,iter)),0,2),3);
                    y=nanmean(zscore(gdata(:,keptT,subjs_g2(:,iter)),0,2),3);
                    
                    r(ni,:,iter)=circularlagcorr(x',y',lags);
                    
                end
            end
        end
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_gg/' seed  ],'r','lags','networks','keptT','subjs_g2','subjs_g1');
        clear r
    end
end



