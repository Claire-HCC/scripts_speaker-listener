%% ttest across subject
% nopeak lag0  ISC diagnal with this method. 
% even pc1 onlycan only explain around 5%. 
% In this case, maybe only the combination of pcs is meaninigfu and pc1 by
% itself is not.

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

for ei=1%:6;%:6;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/networks2networks' '_binSize' num2str(binSize)  '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT2tail_peaks' ],...
            'networks','rzm','lags','rz','keptT');
        [~,~,tn,~,listenerN]=size(rz);
        
        keptSubj=find(squeeze(sum(~isnan(rz(1,1,:,1,:)),3))>0);
        keptT=find(squeeze(sum(~isnan(rz(1,1,:,1,:)),5))>0);
        
        X=permute(rz(:,:,keptT,:,keptSubj),[1 2 4 3 5]);
        X=reshape(X,size(X,1)*size(X,2)*size(X,3),length(keptT)*length(keptSubj));
        keptTS=~isnan(X(1,:));
        [coeff,score,latent,tsquared,explained,mu] = pca(X(:,keptTS),'centered',1);
        coeff_temp=nan([tn listenerN]);
        for pci=2;
            coeff_temp(keptT,keptSubj)=reshape(coeff(:,pci),length(keptT),length(keptSubj));
            coeff=coeff_temp;
            pc=reshape(score(:,pci),6,6,21);

            for sdi=1:6;
                for ni=1:6;
                    [~,peakLagi]=max((pc(sdi,ni,:)),[],3);
                    peakLags(sdi,ni)=(lags(peakLagi));
                    peaks(sdi,ni)=pc(sdi,ni,peakLagi);
                end
            end
            save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/networks2networks' '_binSize' num2str(binSize)  '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_pcs' num2str(pci) ],...
                'pc','peaks','peakLags','networks','keptT','keptSubj','coeff','score','latent','tsquared','explained','mu');
        end
    end
end

