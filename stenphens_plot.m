clear all
close all
set_parameters;
lags=-10:10;
figure;
% for ei=1:4;
%     exp=experiments{ei};
%     subplot(2,2,ei)
%   %  load([expdir exp '\fmri\pattern\regression\tr\roi\mor\SL_g\lag' num2str(min(lags)) '-'  num2str(max(lags)) '_permPhase_stats.mat'],'sig_fdr')
%     
%      load([expdir exp '\fmri\pattern\regression\tr\roi\mor\SL_g\lag' num2str(min(lags)) '-'  num2str(max(lags)) ],'b')
%    
%     plot(lags,nanmean(zscore(b(:,2:end),0,2)),'r')
%     hold on
%     load(['Y:\claire\speaker-listener\' exp '\fmri\pattern\regression\tr\roi\mor\LL_leave1out\lag'  num2str(min(lags)) '-'  num2str(max(lags))  '.mat'])
%     plot(lags,nanmean(nanmean(zscore(b(:,2:end,:),0,2),3)),'b')
%     hold off
%     grid on
%     set(gca,'xtick',[min(lags):4:max(lags)],'xticklabel',[min(lags):4:max(lags)]*1.5)
%     ylabel('zscored beta values');
%     title({exp,'spatial-temporal regression','averaged across all ROIs'})
%     xlabel('speaker precedes<----sec---->Listeners precede')
%     set(gca,'fontsize',14)
%    
% end


figure;
for ei=1:4;
    exp=experiments{ei};
    subplot(2,2,ei)
  % load([expdir exp '\fmri\pattern\regression\tr\roi\mor\SL_g\lag' num2str(min(lags)) '-'  num2str(max(lags)) '_permPhase_stats.mat'],'sig_fdr')
    
     load([expdir exp '\fmri\pattern\lagcorr\tr\roi\mor\SL_g\lag' num2str(min(lags)) '-'  num2str(max(lags)) ],'r')
   
    plot(lags,nanmean(zscore(r,0,2)),'r')
    hold on
    load(['Y:\claire\speaker-listener\' exp '\fmri\pattern\lagcorr\tr\roi\mor\LL_leave1out\lag'  num2str(min(lags)) '-'  num2str(max(lags))  '.mat'])
    plot(lags,nanmean(nanmean(zscore(r,0,2),3)),'b')
    hold off
    grid on
    set(gca,'xtick',[min(lags):4:max(lags)],'xticklabel',[min(lags):4:max(lags)]*1.5)
    ylabel('zscored r values');
    title({exp,'spatial-temporal regression','averaged across all ROIs'})
    xlabel('speaker precedes<----sec---->Listeners precede')
    set(gca,'fontsize',14)
   
end

