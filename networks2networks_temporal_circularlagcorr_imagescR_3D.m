close all
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='restFc_isc30PercMasked_75Overlap_cluster6';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};

lags=-15:15;
exps={'sherlock','merlin','pieman_old','black','forgot','ABC','bronx','pieman','pieman_oldWord','pieman_rest','crossExps'};
eis=[1:11];
fs=[];

for ei=1%:length(eis);
    exp=exps{eis(ei)};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_lag-15-15_timeReversed_peaks' ],...
        'keptT','networks','tmid','rzm','peakLags','pfdr_peaks','peaks','npeaks');
    
    rzm=zscore(rzm,0,3);
    
    surfmat=zeros(size(permute(rzm(:,:,[-15:15]+tmid),[1 3 2]))+1);
    surfmat(1:(end-1),1:(end-1),1:(end-1))=permute(rzm(:,:,[-15:15]+tmid),[1 3 2]);
    
    fsize=[16.5 7];
    
    clim=[1 3.5];
    cols=gray(100);
    cols=cols(discretize(log2([clim(1):0.01:clim(2)]+1),100),:);
    
    figure('unit','centimeter','position',[0 0 fsize]);
    subplot(1,2,1);
    h=slice(surfmat(:,1:15,:),1:15,[],[]);
    set(h,'edgecolor','none');
    axis off
    set(gca,'YDir','reverse','Zdir','reverse')
    colormap(cols);
    caxis(clim)
    
    
    subplot(1,2,2);
    h=slice(surfmat(:,16:31,:),1:15,[],[]);
    set(h,'edgecolor','none');
    axis off
    set(gca,'YDir','reverse','Zdir','reverse')
    colormap(cols);
    caxis(clim)
    
    
    fsize=[ 7 7];
    figure('unit','centimeter','position',[0 0 fsize]);
    
    im=imagesc(squeeze(rzm(:,:,tmid)),clim);
    lm=get(gca,'xlim');
    colormap(cols)
    set(gca,'xtick',1:6,'ytick',1:6,'color','k','xlim',lm,'ylim',lm);
    for ni=1:6;
        pos=[ni-0.5 ni-0.5 1 1];
        rectangle('Position',pos ,'edgecolor','k','linewidth',1.5);
    end
    set(gca,'ytick',1:length(networks),'yticklabels',networks);
    set(gca,'xtick',1:length(networks),'xticklabels',networks);
    xtickangle(75);
    ylabel('Network','fontsize',16,'fontweight','bold')
    xlabel('Network','fontsize',16,'fontweight','bold')
    
    
    fsize=[13 4.5];
    figure('unit','centimeter','position',[0 0 fsize]);
    im=imagesc(squeeze(rzm(6,:,tmid+[-15:15])),clim);
    lm=get(gca,'xlim');
    colormap(cols);
    set(gca,'xtick',[1 16 31],'xticklabels',strrep(cellstr(num2str(1.5*[-15 0 15]')),' ',''));
    ylabel([networks{6} ' seed'],'FontWeight','bold','fontsize',16)
    xlabel('Lag (sec)','FontWeight','bold','fontsize',16,'fontweight','bold');
    set(gca,'ytick',1:length(networks),'yticklabels',networks);
    
    
    figure('unit','centimeter','position',[0 0 fsize]);
    rh = rectangle('position',[ ...
        im.XData(1)-.5,...
        im.YData(1)-.5,...
        range(im.XData)+1,...
        range(im.YData)+1],...
        'EdgeColor', 'k','linewidth',2);
    axis off
end