clear all
close all

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
lags_tested={-10:10, -20:20, -30:30, -10:-4, -20:-4, -30:-4, -10:-1};

figure
subplot_ind=reshape(1:8,4,2)';
for lagi=1%:length(lags_tested);
    lags=lags_tested{lagi};
    
    for ei=1:4;
        exp=experiments{ei};
        
        load([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLeach/regression_SLeach_lag' num2str(min(lags)) '-' num2str(max(lags))  ],'b');
        b_s=b(:,2:end,:);
        subjn=size(b_s,3);
        b_s_m=nanmean(b_s,3);
        
        b_l=nan([size(b_s) subjn]);
        for perm=1:subjn;
            load(sprintf('%s/%s/fmri/temporal_regression/%s/roi/%s/SLeach/perm/regression_SLeach_lag%d-%d_permSL%03d',expdir,exp,timeUnit,froidir,min(lags),max(lags),perm),'b');
            b_l(:,:,perm,:)=reshape(b(:,2:end,:),[size(b,1) length(lags) 1 subjn]);
        end
        
        b_l_m=nanmean(nanmean(b_l,3),4);
        
        ris=find(sum(isnan(b_s_m),2)==0);
        
        subplot(2,4,subplot_ind(1,ei));
        imagesc(zscore(b_l_m(ris,:),0,2));
    %    line([find(lags==0) find(lags==0)],[ 0 sum(sig_fdr)+0.5],'color','w');
        title({exp,['Listener-Listener Weightings(z)']});
        set(gca,'ytick',[]);
        set(gca,'xtick',1:5:length(lags),'xticklabels',lags(1:5:end));
        %  set(gca,'ytick',1:length(ris),'yticklabels',rnames(ris),'fontsize',6);
        if ei==1;
            xlabel('Speaker precedes-----TR(1.5s)-----Listeners precede');
            ylabel(['ROIs']);
        end
        
        
        subplot(2,4,subplot_ind(2,ei));
        imagesc(zscore(b_s_m(ris,:),0,2));
        set(gca,'xtick',1:5:length(lags),'xticklabels',lags(1:5:end));
        %   set(gca,'ytick',1:length(ris),'yticklabels',rnames(ris),'fontsize',6);
        
        % line([find(lags==0) find(lags==0)],[ 0 sum(sig_fdr)+0.5],'color','w')
        title(['Speaker-Listener Weightings(z)'])
        
        clear b_l b_s
    end
    
end

%             figure
%             boxplot(acc')
%             set(gca,'xtick',1:length(rnames),'xticklabels',rnames);
%             xtickangle(75)
%             grid on

% %
% figure; ciplot_claire(zscore(squeeze(b_s(ri,:,:))',0,2),lags,'r',0.3);
% hold on;
% ciplot_claire(zscore(squeeze(b_l(ri,:,:))',0,2),lags,'b',0.3)
% title(rnames(ri))
