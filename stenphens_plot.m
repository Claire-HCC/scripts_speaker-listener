binSize=10;
exp='sherlock'
load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
r2_sl=r2;
load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/binSize' num2str(binSize) '_lag0-0' ],'r2','rnames');
r2_ll=r2;
r2_llm=nanmean(r2_ll,3);
figure; scatter(r2_sl(1,:),r_llm(:),40,'k','fileld')
figure; scatter(r2_sl(1,:),r2_llm(:),40,'k','fileld')
figure; scatter(r2_sl(1,:),r2_llm(1,:),40,'k','fileld')
figure; subplot(2,2,1); scatter(r2_sl(1,:),r2_llm(1,:),40,'k','fileld')
subplot(2,2,1); plot(r2_sl(1,:),r2_llm(1,:),40,'k','fileld')
subplot(2,2,1); plot(r2_sl(1,:),r2_llm(1,:))
binSize=40;
load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
r2_sl=r2;
load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/binSize' num2str(binSize) '_lag0-0' ],'r2','rnames');
r2_ll=r2;
r2_llm=nanmean(r2_ll,3);
subplot(2,2,3); plot(r2_sl(1,:),r2_llm(1,:),40,'k','fileld')
subplot(2,2,3); scatter(r2_sl(1,:),r2_llm(1,:),40,'k','fileld')
subplot(2,2,1); plot(r2_sl(1,:),r2_llm(1,:))
subplot(2,2,4); plot(r2_sl(1,:),r2_llm(1,:))
binSize=10;
subplot(2,2,1); plot(r2_sl(1,:),r2_llm(1,:),40,'k','fileld')
load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
r2_sl=r2;
load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/binSize' num2str(binSize) '_lag0-0' ],'r2','rnames');
r2_ll=r2;
r2_llm=nanmean(r2_ll,3);

subplot(2,1,1);
subplot(2,1,2); plot([r2_sl(1,:); r2_llm(1,:)]')

% 
% 
% fsize=[10 10]
% for ei=3:4;
%     exp=experiments{ei};
%     
%    
%    load(['Y:\claire\speaker-listener\' exp '\fmri\pattern\lagcorr\tr\roi\mor\LL_leave1out\lag'  num2str(min(lags)) '-'  num2str(max(lags))  '.mat'])
%     figure('unit','centimeter','position',[0 0 fsize]);
% 
%      ciplot_claire(nanmean(r,3),lags,'r',0.3);
%   %  plot(lags,nanmean(atanh(r),3),'k')
%     grid on
%     set(gca,'xtick',[min(lags):5:max(lags)],'xticklabel',[min(lags):5:max(lags)])
%     ylabel('R(z)');
%     title(exp)
%     xlabel('Lag (TR=1.5 sec)')
%     set(gca,'fontsize',14)
%     
% hold on
%      load([expdir exp '\fmri\pattern\lagcorr\tr\roi\mor\SL_each\lag' num2str(min(lags)) '-'  num2str(max(lags)) ],'r')
% 
%      ciplot_claire(nanmean(r,3),lags,'k',0.3);
%      hold off
%    % plot(lags,nanmean(atanh(r),3),'k')
%     grid on
%     set(gca,'xtick',[min(lags):5:max(lags)],'xticklabel',[min(lags):5:max(lags)])
%     ylabel('R(z)');
%     title(exp)
%     xlabel('Speaker precedes---Lag(TR)---Listener precedes','fontsize',12)
%     set(gca,'fontsize',14)
% end
% 
% 
