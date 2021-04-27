clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rgb=roi_table.rgb;
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags=-20:20
binSize=30;

network_newOrder={
    'Auditory_Language',...
    'DMN2',...
    'Attention',...
    'Executive',...
    'DMN1',...
    'Visual'...
    };
cols=jet(7);
gray=[0.9 0.9 0.9];
eis=[1 2 4 11 12  13];%
fsize=[45 7];
figure('unit','centimeter','position',[0 0 fsize],'color','k');

for i=1:length(eis);
    ei=eis(i);
    exp=experiments{ei};
    subplot(1,6,i)
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/rois2rois_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaksMax' ],...
        'rnames','rzm','rz','lags','keptT','pfwe','peaks_pfwe','peakLags_pfwe','pfdr','peaks_p','peakLags_p','peakLags','peaks');
    
    for ri=1:length(rnames);
        
        line([min(lags) max(lags)],[0 0],'color',gray)
        hold on;
        
        temp=squeeze(rz(ri,ri,:,:))';
        
        network=roi_table.network(ri);
        ni=find(ismember(network_newOrder,network));
      %  ciplot_claire(temp,lags*tr(ei),[cols(ni+1,:) ],0.1);
        plot(lags*tr(ei),nanmean(temp)/max(nanmean(temp)),'color',cols(ni+1,:))
        hold on;
        
        set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)
        %     ylim([-0.2 0.45])
        %  xlim([min(lags) max(lags)]);
        xlim([-3 3])
        
        title(strrep(exp,'_',' '),'fontsize',14,'color','w');
        % set(gca,'xtick',[-20:1:20],'xticklabel',[-20:1:20]*tr(ei));
        xlabel('Lag (sec)');
        ylabel('Normalized R (z)')
        set(gca,'fontsize',10)
        % grid on
        
    end
    
end
