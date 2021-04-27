%% ttest across subject
% fdr roix lag or fdr roixroixlag?
close all
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
lags_tested={-10:10,  -40:40};
binSize=30;

for ei=1:6;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/networks2networks' '_binSize' num2str(binSize)  '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT2tail_peaks' ],...
            'networks','rzm','lags','rz','keptT');
        [~,~,tn,~,listenerN]=size(rz);
        
        for pci=1;
            coeffs=nan([tn listenerN]);
            pcs=nan([6 6 21 listenerN]);
            explaineds=nan([listenerN 1]);
            for si=1:listenerN;
                if sum(~isnan(squeeze(rz(1,1,:,1,si))))>0;
                    X=permute(rz(:,:,:,:,si),[1 2 4 3]);
                    X=reshape(X,size(X,1)*size(X,2)*size(X,3),size(X,4));
                    keptT=~isnan(X(1,:));
                    [coeff,score,latent,tsquared,explained,mu] = pca(X(:,keptT));
                    coeffs(keptT,si)=coeff(:,pci);
                    pcs(:,:,:,si)=reshape(score(:,pci),6,6,21);
                    explaineds(si)=explained(pci);
                else;
                    coeffs(:,si)=NaN;
                    pcs(:,:,:,si)=NaN;
                    explaineds(pci,si)=NaN;
                end
            end
            
            pcm=nanmean(pcs,4);
            peakLags=nan(6,6);
            peaks=nan(6,6);
            for sdi=1:6;
                for ni=1:6;
                    [~,peakLagi]=max((pcm(sdi,ni,:)),[],3);
                    peakLags(sdi,ni)=(lags(peakLagi));
                    peaks(sdi,ni)=pcm(sdi,ni,peakLagi);
                    [~,p(sdi,ni,:),~,stats]=ttest(squeeze(pcs(sdi,ni,:,:))');
                    t(sdi,ni,:)=stats.tstat;
                end
            end
            save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/networks2networks' '_binSize' num2str(binSize)  '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_pcs' num2str(pci) ],...
                'pcs','peaks','peakLags','networks','coeffs','explaineds','keptT');
            %             [~,~,pfdr]=fdr(p(:));
            %             pfdr=reshape(pfdr,size(p));
            %             pcm_temp=pcm;
            %             pcm_temp(pfdr>.05)=NaN;
            %             peakLags_pfdr=nan(6,6);
            %             peaks_pfdr=nan(6,6);
            %
            %             for sdi=1:6;
            %                 for ni=1:6;
            %                     if nansum(~isnan(squeeze(squeeze(pcm_temp(sdi,ni,:)))))>0;
            %                         [~,peakLagi]=max((pcm_temp(sdi,ni,:)),[],3);
            %                         peakLags_pfdr(sdi,ni)=(lags(peakLagi));
            %                         peaks_pfdr(sdi,ni)=pcm_temp(sdi,ni,peakLagi);
            %                     end
            %                 end
            %             end
            %
            %             save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/networks2networks' '_binSize' num2str(binSize)  '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_pcs' num2str(pci) ],...
            %                 'pcs','peaks','peakLags','networks','peakLags_pfdr','peaks_pfdr','coeffs');
        end
    end
end

