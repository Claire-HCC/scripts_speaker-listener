%% ttest across subject
% fdr roix lag or fdr roixroixlag?

clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
lags_tested={-20:20,  -40:40};
binSize_tested=[ 30 20 100];

for ei=9:10;%[1 2 4  11:12];
    exp=experiments{ei};
    
    for bi=[1]%:2;
        binSize=binSize_tested(bi);
        
        for lagi=1%:length(lags_tested);
            lags=lags_tested{lagi};
            
            load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/' networks{1} '_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'networks','r','lags','keptT');
            [~,tn,~,listenerN]=size(r);
            rz=nan([length(networks) length(networks) tn length(lags) listenerN]);
            
            for sdi=1:length(networks);
                seed=networks{sdi};
                if exist([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/' seed '_binSize' num2str(binSize)  '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat']);
                    load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/' seed  '_binSize' num2str(binSize)  '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'networks','r','lags','keptT');
                    rz(sdi,:,:,:,:)=atanh(r);
                end
            end
            
            rzm=nanmean(rz,5);
            p=nan(size(rzm));
            t=p;
            
            for sdi=1:length(networks);
                seed=networks{sdi};
                for ni=1:length(networks);
                    for ti=1:tn;
                        if ~isnan(r(1,ti,1,1));
                            [~,p(sdi,ni,ti,:),~,stats]=ttest(squeeze(rz(sdi,ni,ti,:,:))',0,'tail','both');
                            t(sdi,ni,ti,:)=stats.tstat;
                        end
                    end
                end
            end
            
            pfwe=p*(sum(~isnan(p(:))));
            
            pfdr=nan(size(p(:)));
            [ ~,~,pfdr(~isnan(p(:)))]=fdr(p(~isnan(p(:))));
            pfdr=reshape(pfdr,size(p));
            
            peaks=nan([length(networks) length(networks) tn]);
            peakLags=peaks;
            peaks_pfwe=nan([length(networks) length(networks) tn]);
            peakLags_pfwe=peaks_pfwe;
            peaks_pfdr=nan([length(networks) length(networks) tn]);
            peakLags_pfdr=peaks_pfdr;
            rz_tempfwe=rzm;
            rz_tempfwe(pfwe>.05)=NaN;
            
            rz_tempfdr=rzm;
            rz_tempfdr(pfdr>.05)=NaN;
            
            for sdi=1:length(networks);
                for ni=1:length(networks);
                    for ti=1:tn;
                        
                        if sum(isnan(rzm(sdi,ni,ti,:))==0);
                            
                            [~, peakLagi]=max((rzm(sdi,ni,ti,:)),[],4);%max(abs(rzm(ni,:)),[],2);
                            peakLags(sdi,ni,ti)=(lags(peakLagi));
                            peaks(sdi,ni,ti)=rzm(sdi,ni,ti,peakLagi);
                            
                            if min(pfdr(sdi,ni,:))<.05;
                                [~, peakLagi]=max((rz_tempfdr(sdi,ni,ti,:)),[],4);;%max(abs(rz_tempfdr(sdi,ni,:)),[],2);
                                peakLags_pfdr(sdi,ni,ti)=(lags(peakLagi));
                                peaks_pfdr(sdi,ni,ti)=rz_tempfdr(sdi,ni,ti,peakLagi);
                            end
                            
                            if min(pfwe(sdi,ni,:))<.05;
                                [~, peakLagi]=max((rz_tempfwe(sdi,ni,ti,:)),[],4);%max(abs(rz_tempfwe(sdi,ni,:)),[],2);
                                peakLags_pfwe(sdi,ni,ti)=(lags(peakLagi));
                                peaks_pfwe(sdi,ni,ti)=rz_tempfwe(sdi,ni,ti,peakLagi);
                            end
                        end
                    end
                end
            end
            save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/networks2networks' '_binSize' num2str(binSize)  '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT2tail_peaks' ],...
                'networks','t','rz','rzm','lags','keptT','p','pfwe','peaks_pfwe','peakLags_pfwe','pfdr','peaks_pfdr','peakLags_pfdr','peakLags','peaks');
        end
    end
end

