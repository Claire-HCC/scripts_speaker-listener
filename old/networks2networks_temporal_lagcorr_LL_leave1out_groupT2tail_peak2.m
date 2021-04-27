%% find the peak nearest to lag 0 instead of the absolute peak
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
lags_tested={-20:20 , -15:15, -10:10,  -40:40};

for ei=[1 2 4 11 12 9 10];
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks_peakLag0_mask/' froidir '/LL_leave1out/perm/nwtworks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'networks','r','lags','keptT');
        [~,~,listenerN,~,~]=size(r);
        rzm_perm=nanmean(rz,4);
        
        for sdi=1:length(networks);
            seed=networks{sdi};
            if exist([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks_peakLag0_mask/' froidir '/LL_leave1out/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat']);
                load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks_peakLag0_mask/' froidir '/LL_leave1out/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'networks','r','lags','keptT');
                rz(sdi,:,:,:)=atanh(r);
            end
        end
        rzm=nanmean(rz,4);
        
        p=mean(rzm<rzm_perm,3);
pfwe=p*(size(rzm,1)*size(rzm,2)*size(rzm,3));
% plot(lags,rzm');
% for ni=1:6;
% text(lags(pfwe(ni,:)<.05),rzm(ni,pfwe(ni,:)<.05),'*');
% end



        peakLags=nan([length(networks), length(networks) ]);
        peaks=nan([length(networks), length(networks) ]);
        peaks_pfdr=nan([length(networks) length(networks)]);
        peakLags_pfdr=peaks_pfdr;
        rzm_temp=rzm;
        % only consider positve R
        rzm_temp(rzm_temp<0)=NaN;
        
        rzm_temppfdr=rzm_temp;
        rz_tempfdr(pfdr>.05)=NaN;
        for sdi=1:6;
            for ni=1:6;
                [pks, locs]=findpeaks(squeeze(rzm_temp(sdi,ni,:)));
                if ~isempty(pks);
                    [~,loci]=(min(abs(locs-find(lags==0))));
                    peakLags(sdi,ni)=lags(locs(loci));
                    peaks(sdi,ni)=pks(loci);
                end
                
                [pks, locs]=findpeaks(squeeze(rzm_temppfdr(sdi,ni,:)));
                if min(pfdr(sdi,ni,:))<.05;
                    [~,loci]=(min(abs(locs-find(lags==0))));
                    peakLags_pfdr(sdi,ni)=lags(locs(loci));
                    peaks_pfdr(sdi,ni)=pks(loci);
                end
            end
        end
        
        peakLags_subj=nan([length(networks), length(networks) listenerN ]);
        peaks_subj=nan([length(networks), length(networks) listenerN]);
        rz_temp=rz;
        % only consider positve R
        rz_temp(rz_temp<0)=NaN;
        for sdi=1:6;
            for ni=1:6;
                for si=1:listenerN;
                    [pks, locs]=findpeaks(squeeze(rz_temp(sdi,ni,:,si)));
                    if ~isempty(pks);
                        [~,loci]=(min(abs(locs-find(lags==0))));
                        peakLags_subj(sdi,ni,si)=lags(locs(loci));
                        peaks_subj(sdi,ni,si)=pks(loci);
                    end
                    
                    
                end
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks_peakLag0_mask/' froidir '/LL_leave1out/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT2tail_peaksNearest2zero' ],...
            'networks','t','rz','rzm','lags','keptT','p','pfwe','pfdr','peakLags','peaks','peakLags_subj','peaks_subj','peaks_pfdr','peakLags_pfdr');
    end
end

