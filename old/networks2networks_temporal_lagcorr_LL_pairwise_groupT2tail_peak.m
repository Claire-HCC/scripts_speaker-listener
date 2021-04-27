%% ttest across subject
% fdr roix lag or fdr roixroixlag?

clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
lags_tested={-10:10,  -40:40};

for ei=9;%[3 5 7 8 9:12];
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_pairwise/' networks{1} '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'networks','r','lags','keptT');
        [~,~,listenerN,~]=size(r);
        rz=nan([length(networks) length(networks) length(lags) listenerN listenerN]);
        
        for sdi=1:length(networks);
            seed=networks{sdi};
            if exist([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_pairwise/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat']);
                load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_pairwise/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'networks','r','lags','keptT');
                rz(sdi,:,:,:,:)=atanh(r);
            end
        end
        rz=reshape(rz,length(networks), length(networks) ,length(lags), listenerN*listenerN);
        rz_temp=rz;
        rz_temp(rz_temp<0)=NaN;
        peakLags=nan([length(networks), length(networks) , size(rz,[ 4])]);
        peaks=nan([length(networks), length(networks) , size(rz,[ 4])]);
        for sdi=1:length(networks);
            for ni=1:length(networks);
                for spi=1:size(rz,4);
                  
                    [pks, locs]=findpeaks(squeeze(rz(sdi,ni,:,spi)),'threshold',0);
                    if ~isempty(pks);
                    [~,loci]=(min(abs(locs-find(lags==0))));
                    peakLags(sdi,ni,spi)=lags(locs(loci));
                    peaks(sdi,ni,spi)=pks(loci);
                    end
                end
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_pairwise/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT2tail_peaks' ],...
            'networks','rz','lags','keptT','peakLags','peaks');
    end
end

