clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';
binSize=[10]; % tr;
lags=-10:10  ;

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

ri=1;
rname=rnames{ri};
for ei=1:4;
    exp=experiments{ei};
    load([ expdir exp '\fmri\pattern\regression\tr\roi\mor\SL_g_bined\binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'b')
fsize=[22 9];
figure('unit','centimeter','position',[0 0 fsize]);
subplot(1,2,1);
imagesc(squeeze(b(ri,:,2:end)),[-0.08 0.08])
colorbar
title({[upper(exp(1)) exp(2:end) ', ' rname],['binSize' num2str(binSize) ', lag' num2str(min(lags)) '~' num2str(max(lags)) ]})
xlabel('Lag (TR)');
ylabel('Time (TR)');
set(gca,'fontsize',14,'xtick',1:5:length(lags),'xticklabels',lags(1:5:length(lags)));

subplot(1,2,2);
ciplot_claire(squeeze(b(ri,:,2:end)),lags,'k',0.3);
title({[upper(exp(1)) exp(2:end) ', ' rname],['binSize' num2str(binSize) ', lag' num2str(min(lags)) '~' num2str(max(lags)) ]});
grid on
xlabel('Lag (TR)');
ylabel('beta-value');
set(gca,'fontsize',14);

end