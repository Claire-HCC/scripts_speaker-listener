%% find the peak nearest to lag 0 instead of the absolute peak
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
lags_tested={-20:20 , -15:15, -10:10,  -40:40};

for ei=[1 2 4 9:12];
    exp=experiments{ei};
    
    for lagi=3%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks_peakLag0_mask/' froidir '/LL_leave1out/perm/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'networks','r','lags','keptT');
        [~,~,listenerN,~,~]=size(r);
        rzm_perm=max((squeeze(nanmean(atanh(r),4))),[],3);
        rzm_perm=repmat(rzm_perm,1,1,length(lags),1);
        
        rz=[];
        for sdi=1:length(networks);
            seed=networks{sdi};
            if exist([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks_peakLag0_mask/' froidir '/LL_leave1out/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat']);
                load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks_peakLag0_mask/' froidir '/LL_leave1out/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'networks','r','lags','keptT');
                rz(sdi,:,:,:)=atanh(r);
            end
        end
        rzm=nanmean(rz,4);
        
        p=mean(rzm<rzm_perm,4);
        pfwe=p*(size(rzm,1)*size(rzm,2));
        
        peakLags=nan([length(networks), length(networks) ]);
        peaks=nan([length(networks), length(networks) ]);
        peaks_pfwe=nan([length(networks) length(networks)]);
        peakLags_pfwe=peaks_pfwe;
        rzm_temp=rzm;
        % only consider positve R
        rzm_temp(rzm_temp<0)=NaN;
        
        rzm_temppfwe=rzm_temp;
        rz_tempfwe(pfwe>.01)=NaN;
        for sdi=1:6;
            for ni=1:6;
                [pks, locs]=findpeaks(squeeze(rzm_temp(sdi,ni,:)));
                if ~isempty(pks);
                    %    [~,loci]=(min(abs(locs-find(lags==0))));
                    [~,loci]=max(pks);
                    peakLags(sdi,ni)=lags(locs(loci));
                    peaks(sdi,ni)=pks(loci);
                end
                
                [pks, locs]=findpeaks(squeeze(rzm_temppfwe(sdi,ni,:)));
                if ~isempty(pks) & min(pfwe(sdi,ni,:))<.05;
                    %    [~,loci]=(min(abs(locs-find(lags==0))));
                    [~,loci]=max(pks);
                    peakLags_pfwe(sdi,ni)=lags(locs(loci));
                    peaks_pfwe(sdi,ni)=pks(loci);
                end
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks_peakLag0_mask/' froidir '/LL_leave1out/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_peaks' ],...
            'networks','rzm_perm','rzm','lags','keptT','p','pfwe','peakLags','peaks','peaks_pfwe','peakLags_pfwe');
    end
end

