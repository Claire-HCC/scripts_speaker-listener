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
lags_tested={-20:20,  -40:40};
binSize=30;
% network_newOrder={
%     'Auditory_Language',...
%     'DMN2',...
%     'Attention',...
%     'Executive',...
%     'DMN1',...
%     'Visual'...
%     };

for ei=[1 2 4 11 12 9 10];
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/networks2networks' '_binSize' num2str(binSize)  '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT2tail_peaks' ],...
            'networks','rzm','lags');
        
        X=permute(rzm,[1 2 4 3]);
        X=reshape(X,size(X,1)*size(X,2)*size(X,3),size(X,4));
        keptT=~isnan(X(1,:));
        [coeff,score,latent,tsquared,explained,mu] = pca(X(:,keptT),'centered','on');
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/networks2networks' '_binSize' num2str(binSize)  '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_pc'  ],...
            'networks','coeff','score','latent','tsquared','explained','mu');
        fsize=[28 8.5];
        figure('unit','centimeter','position',[0 0 fsize]);
        subplot(1,3,1);
        plot(find(keptT),coeff(:,1:2),'linewidth',2);
        legend('pc1','pc2'); legend boxoff
        set(gca,'fontsize',12);
        grid on;
        set(gca,'xticklabels',get(gca,'xtick')*tr(ei))
        xlabel('Time (sec)')
        ylabel('Component Weighting');
        title({[upper(exp(1)) exp(2:end)],sprintf('Explained variance: %d%%+%d%%',round(explained(1)),round(explained(2)))});
        
        for pci=1:2;
            pc=reshape(score(:,pci),6,6,length(lags));
            pc_temp=pc;
            pc_temp(pc_temp<0)=NaN;
            peakLags=nan([length(networks), length(networks) ]);
            peaks=nan([length(networks), length(networks) ]);
            for sdi=1:6;
                for ni=1:6;
                    [pks, locs]=findpeaks(squeeze(pc_temp(sdi,ni,:)));
                    if ~isempty(pks);
                        [~,loci]=(min(abs(locs-find(lags==0))));
                        peakLags(sdi,ni)=lags(locs(loci));
                        peaks(sdi,ni)=pks(loci);
                    end
                end
            end
            save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/networks2networks' '_binSize' num2str(binSize)  '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_pc' num2str(pci) '_peaksNearest2zero' ],...
                'pc','peaks','peakLags','networks');
            
            [~,nis]=sort(nanmean(peakLags));
            network_newOrder=networks(nis);
            % nis=cellfun(@(x) find(ismember(networks,x)),network_newOrder);
            subplot(1,3,pci+1);
         
               peakLags_bined= discretize(peakLags,   [-20 -15 -10 -6 -3 0 1 4  7  11 16 20]);
            imagesc(peakLags_bined(nis,nis),[1 11]);
          %  peakLags_bined= discretize(peakLags,[-10 -6 -3 0 1 4  7  10]);
           % imagesc(peakLags_bined(nis,nis),[1 7]);
            colormap jet
            title(['pc' num2str(pci)])
            for ni=1:length(network_newOrder);
                network=network_newOrder{ni};
                temp=find(ismember(network_newOrder,network));
                
                pos=[min(temp)-0.5 min(temp)-0.5 length(temp) length(temp)];
                rectangle('Position',pos ,'edgecolor','w','linewidth',0.5);
                
                if strcmp(network,'Auditory_Language'); network={'Auditory/','Language'};end
                text(pos(1)+pos(3)/2,pos(2)+pos(4)*0.3,network,'HorizontalAlignment','center','color','w','fontweight','bold','fontsize',9);
                
            end
        end
    end
end
pci=1;
ind=triu(ones(6,6),1);
clear pl; for ei=1:8;exp=experiments{ei};
    load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/networks2networks' '_binSize' num2str(binSize)  '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_pc' num2str(pci) ],...
        'pc','peaks','peakLags');temp=peakLags(ind==1);
    pl(:,ei)=temp(:);
end
[r p]=corr(pl,'type','spearman')
