function network2networks_temporal_lagcorr_LL_leave1out_perm(i)
eis=[1 2 4 9:12];
ei=eis(i);
loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
seeds=networks;
crop_start=25;
crop_end=20;
lags_tested={-20:20, -15:15, -10:10, -40:40, -60:60};
permN=1000;
exp=experiments{ei};

mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks_peakLag0_mask/' froidir '/LL_leave1out/perm/']);
for lagi=3;%;%:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' networks{1} '_peakLag0_mask.mat' ],'gdata');
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    
    gdatas=[];
    for ni=1:size(networks);
        network=networks{ni};
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '_peakLag0_mask.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '_peakLag0_mask.mat' ],'gdata');
            
            gdata(:,:,subjects_excluded{ei})=NaN;
            gdata=zscore(gdata(:,keptT,:),0,2);
            gdata=nanmean(gdata,1);
            gdatas(ni,:,:)=gdata;
        end
    end
    
    r=nan([length(networks)  length(networks)  length(lags) listenerN  permN]);
    for perm=1:permN;
        rng(perm)
        
        gdata_seed_perm=[];
        gdatas_perm=[];
        for si=1:listenerN;
            gdatas_seed_perm(:,:,si)=(phase_rand2(gdatas(:,:,si)',1))';
            gdatas_perm(:,:,si)=(phase_rand2(gdatas(:,:,si)',1))';
        end
        
        for si=1:listenerN;
            othersi=1:listenerN;
            othersi=othersi(othersi~=si);
            
            y=nanmean(zscore(gdatas_perm(:,:,othersi),0,2),3)';
            x=zscore(gdatas_seed_perm(:,:,si),0,2)';
            
            for sdi=1:size(x,2);
                x_temp=repmat(x(:,sdi),1,length(networks));
                
                r(sdi,:,:,si,perm)=lagcorr_claire(x_temp,y,lags)';
            end
        end
        disp(perm)
    end
    
    save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks_peakLag0_mask/' froidir '/LL_leave1out/perm/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags))  ],'r','lags','networks','keptT','-v7.3');
    clear r
end



