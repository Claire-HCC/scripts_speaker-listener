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

for ei=1:4;%1:4;%1:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/SL_each/' networks{1} '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'networks','r','lags','keptT');
        [~,~,listenerN]=size(r);
        rz=nan([length(networks) length(networks) length(lags) listenerN]);
        
        for sdi=1:length(networks);
            seed=networks{sdi};
            if exist([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/SL_each/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat']);
                load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/SL_each/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'networks','r','lags','keptT');
                rz(sdi,:,:,:)=atanh(r);
            end
        end
        
        rzm=nanmean(rz,4);
        p=nan(size(rzm));
        t=p;
        
        for sdi=1:length(networks);
            seed=networks{sdi};
            for ni=1:length(networks);
                [~,p(sdi,ni,:),~,stats]=ttest(squeeze(rz(sdi,ni,:,:))',0,'tail','both');
                t(sdi,ni,:)=stats.tstat;
            end
        end
        
        pfwe=p*(sum(~isnan(p(:))));
        
        pfdr=nan(size(p(:)));
        [ ~,~,pfdr(~isnan(p(:)))]=fdr(p(~isnan(p(:))));
        pfdr=reshape(pfdr,size(p));
        
        peaks=nan([length(networks) length(networks)]);
        peakLags=peaks;
        peaks_pfwe=nan([length(networks) length(networks)]);
        peakLags_pfwe=peaks_pfwe;
        peaks_pfdr=nan([length(networks) length(networks)]);
        peakLags_pfdr=peaks_pfdr;
        rz_tempfwe=rzm;
        rz_tempfwe(pfwe>.05)=NaN;
        
        rz_tempfdr=rzm;
        rz_tempfdr(pfdr>.05)=NaN;
        
        for sdi=1:length(networks);
            for ni=1:length(networks);
                if sum(isnan(rzm(sdi,ni,:))==0);
                    
                    [~, peakLagi]=max(abs(rzm(sdi,ni,:)),[],3);%max(abs(rzm(ni,:)),[],2);
                    peakLags(sdi,ni)=(lags(peakLagi));
                    peaks(sdi,ni)=rzm(sdi,ni,peakLagi);
                    
                    if min(pfdr(sdi,ni,:))<.05;
                        [~, peakLagi]=max(abs(rz_tempfdr(sdi,ni,:)),[],3);;%max(abs(rz_tempfdr(sdi,ni,:)),[],2);
                        peakLags_pfdr(sdi,ni)=(lags(peakLagi));
                        peaks_pfdr(sdi,ni)=rz_tempfdr(sdi,ni,peakLagi);
                    end
                    
                    if min(pfwe(sdi,ni,:))<.05;
                        [~, peakLagi]=max(abs(rz_tempfwe(sdi,ni,:)),[],3);%max(abs(rz_tempfwe(sdi,ni,:)),[],2);
                        peakLags_pfwe(sdi,ni,1)=(lags(peakLagi));
                        peaks_pfwe(sdi,ni,1)=rz_tempfwe(sdi,ni,peakLagi);
                    end
                end
            end
        end
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/SL_each/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT2tail_peaks' ],...
            'networks','t','r','rzm','lags','keptT','p','pfwe','peaks_pfwe','peakLags_pfwe','pfdr','peaks_pfdr','peakLags_pfdr','peakLags','peaks');
    end
end

