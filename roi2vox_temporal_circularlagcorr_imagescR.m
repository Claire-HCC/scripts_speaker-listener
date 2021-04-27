fsize=[16 11];
figure('unit','centimeter','position',[0 0 fsize])
[~,vis]=sort(peakLags);
imagesc(zscore(rzmsm(vis(sig(vis)==1),:),0,2),[0 5.2]); xlim([tmid+[-15 15]]);
colormap hot
set(gca,'xtick',[fliplr(tmid-(15:15:tmid)) tmid:15:tn]);
set(gca,'xticklabels',(get(gca,'xtick')-tmid)*exp_parameters.tr(ei));
xlabel('Lag (sec)')
ylabel('Voxels ordered by lag')
set(gca,'ytick',[])
set(gca,'fontsize',16)
c=colorbar;
c.Label.String = {'Lag-R (z)','normalized across lags'};
title('Precuneus Seed','fontsize',16)
