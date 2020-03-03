
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');

tic % 15 min
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
lags=[-20:20];
 
fsize=[55 35];
f=figure('unit','centimeter','papersize',[fsize],'position',[0 0 fsize],'paperposition',[0 0 fsize]); 
for ei=1:4;
    exp=experiments{ei};
    rnames=table2array(roi_table(:,3));
    ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
    rnames=rnames(ris);
    
load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_betaStats.mat'],'b_real','b_null','b_p','b_t','lags','rnames','keptT','b_sig_fdr_s','b_sig_fdr_l');

subplot(2,4,ei)
% normalization over lags, so to see consistent patters across rois, even
% though the absolute weightings might change over rois.
imagesc(zscore(b_real,0,2),[-3 3])
% imagesc(b_real,[-0.01 0.01]);
title({exp,'Speaker-Listener weightings'})
set(gca,'xtick',1:5:length(lags),'xticklabels',lags(1:5:end));
set(gca,'ytick',1:length(rnames),'yticklabels',strrep(rnames,'_',' '));
ax=gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 5;
% hold on
% p_boundary=(bwboundaries(b_sig_fdr_s));
% visboundaries(p_boundary,'color','w','linewidth',1);
% hold off

subplot(2,4,4+ei)
imagesc(mean(zscore(b_null,0,2),3),[-3 3])
% imagesc(mean(b_null,3),[-0.005 0.005]);
title('Null weightings')
set(gca,'xtick',1:5:length(lags),'xticklabels',lags(1:5:end));
set(gca,'ytick',1:length(rnames),'yticklabels',strrep(rnames,'_',' '));
ax=gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 5;
% hold on
% p_boundary=(bwboundaries(b_sig_fdr_l));
% visboundaries(p_boundary,'color','w','linewidth',0.05);
% hold off

end
