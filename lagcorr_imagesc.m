close all
set_parameters;
lags=-40:40;
froidir='mor';
timeUnit='tr';

fsize=[45,15];
figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize]);
for ei=1:4;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permSL_peaks'],'pfdr','r','lags','rnames','keptT');
 
%     load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'sig_fdr');
%     regression_mask=sig_fdr;
     
    subplot(2,4,ei);
    %  plot(lags,r,'k','linewidth',2)
    imagesc(zscore(r,0,2),[-3.5 3.5])
    set(gca,'xtick',1:10:length(lags),'xticklabels',min(lags):10:max(lags),'ytick',[]);
    
    ylabel('ROIs')
    title(exp)
    c=colorbar('Location','eastoutside');
    c.Label.String = 'zscored-R';
    set(gca,'fontsize',14);
    
    subplot(2,4,ei+4);
    %  plot(lags,r,'k','linewidth',2)
    r_temp=zscore(r,0,2);
     r_temp(pfdr>.05)=NaN;
    % r_temp(regression_mask==0,:)=NaN;
    imagesc(r_temp,[-3.5 3.5])
    set(gca,'xtick',1:10:length(lags),'xticklabels',min(lags):10:max(lags),'ytick',[]);
    
    ylabel('ROIs')
    title(exp)
    c=colorbar('Location','eastoutside');
    c.Label.String = 'zscored-R';
    set(gca,'fontsize',14);
    
end
