function network2networks_temporal_circularlagregress_L_g

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
crop_start=25;
crop_end=20;
lags=-20:-1;
for ei=[4];%[4 11 9 10];
    exp=experiments{ei};
    
    mkdir([expdir '/' exp '/fmri/temporal/circularlagregress/' timeUnit '/network2networks/' froidir '/L_g/']);
    
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
   
    sig=nan([length(networks) length(networks) length(lags)]);
    
    for tgi=1:length(networks);
        target=networks{tgi};
        
        seedsi=1:size(gdata_orig,1);
        seedsi=seedsi(~ismember(seedsi,tgi));
        
        xx=zscore(nanmean(zscore(gdata_orig(seedsi,keptT,:),0,2),3),0,2);
        x=zscore(nanmean(zscore(gdata_orig(tgi,keptT,:),0,2),3),0,2);
        
        sig(seedsi,tgi,:)=circularlagregress2(x',xx',lags);
        
   
 %     save([expdir '/' exp '/fmri/temporal/circularlagregress/' timeUnit '/network2networks/' froidir '/L_g/sig' ],'b','networks','keptT','lags');
    end
    
end





