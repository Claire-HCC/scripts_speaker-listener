close all
clear all
set_parameters;
lags=-10:10;
froidir='mor';
timeUnit='tr';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table','roi_newOrder');

fsize=[25,9];
figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize]);
for ei=3:4;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/LL_leave1out/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r');
    r=nanmean(atanh(r),3);
    
   % load([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r');
   % r=atanh(r);
    
    subplot(1,2,ei-2);
    imagesc(r,[-0.2 0.2]);
    set(gca,'xtick',1:5:length(lags),'xticklabels',min(lags):5:max(lags),'ytick',[]);
    
    ylabel('ROIs')
    title(exp)
    c=colorbar('Location','eastoutside');
    c.Label.String = 'R(z)';
    set(gca,'fontsize',14);
    xlabel('Lag (TR = 1.5 s)')
   % xlabel('Speaker precedes---Lag (TR)---Listener precedes','fontsize',12)
   
    
    %     subplot(2,4,ei+4);
    %
    %     r_temp=r;
    %     r_temp(pfwe>.025)=NaN;
    %     r_temp(peaks_pfwe<0,:)=NaN;
    %     r_temp(r_temp<0)=NaN;
    %     imagesc(r_temp,[-0.15 0.15])
    %     set(gca,'xtick',1:10:length(lags),'xticklabels',min(lags):10:max(lags),'ytick',[]);
    %
    %     ylabel('ROIs')
    %     title(exp)
    %     c=colorbar('Location','eastoutside');
    %     c.Label.String = 'R(z)';
    %     set(gca,'fontsize',14);
    
end
