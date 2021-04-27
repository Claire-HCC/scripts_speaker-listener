function network2networks_temporal_circularlagregress_LL_leave1out

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
crop_start=25;
crop_end=20;

for ei=[4];%[4 11 9 10];
    exp=experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/temporal/circularlagregress/' timeUnit '/network2networks/' froidir '/LL_leave1out/']);
    
    gdata_orig=[];
    for sdi=1:length(networks);
        seed=networks{sdi};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_iscmasked_' seed '.mat' ],'gdata','keptvox');
        gdata(:,:,subjects_excluded{ei})=NaN;
        gdata=(nanmean(gdata,1));
        gdata_orig(sdi,:,:)=gdata;
    end
    
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
    
    
    for tgi=1:length(networks);
        target=networks{tgi};
        
        b=nan([length(networks)  length(lags) listenerN ]);
        
        for si=1:listenerN;
            othersi=1:listenerN;
            othersi=othersi(othersi~=si);
            
            seedsi=1:size(gdata_orig,1);
             seedsi=seedsi(~ismember(seedsi,tgi));
            
            xx=zscore(nanmean(zscore(gdata_orig(seedsi,keptT,othersi),0,2),3),0,2);
            x=zscore(gdata_orig(tgi,keptT,si),0,2);
            
            b(seedsi,:,si)=circularlagregress(x',xx',lags);
            
        end
        save([expdir '/' exp '/fmri/temporal/circularlagregress/' timeUnit '/network2networks/' froidir '/LL_leave1out/target_' target ],'b','networks','keptT','lags');
    end
    
end





