% close all
clear all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');

lags_tested={-20:20};
roi_newOrder={
    'Auditory_Language',...
    'DMN2',...
    'Attention',...
    'Executive',...
    'DMN1',...
    'Visual'...
    };

eis=[ 1 2 4 11 12 13 ];
fsize=[55 12];
figure('unit','centimeter','position',[0 0 fsize]);
for i=[1:length(eis)];
    ei=eis(i);
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/rois2rois_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaksMax' ],...
            'rois','rz','rzm','lags','keptT','p','pfwe','peaks_pfwe','peakLags_pfwe','pfdr','peaks_p','peakLags_p','peakLags','peaks');
        
        subplot(1,6,i);
        triui=triu(ones(size(peaks)),1);
        ris_sig=find(triui(:)==1  & ~isnan(peaks_pfwe(:)));
        ris_unsig=find(triui(:)==1  & ~isnan(peaks(:)));
        
        scatter(peakLags_pfwe(:)*tr(ei),peaks_pfwe(:),20,'filled','k');
        hold on
        scatter(peakLags(:)*tr(ei),peaks(:),20,'k','filled','MarkerFaceAlpha',0.2);
        
        %    scatter(peaks_pfwe(:),abs(peakLags_pfwe(:)*tr(ei)),20,'filled','k');
        %   hold on
        %  scatter(peaks(:),abs(peakLags(:)*tr(ei)),20,'k','filled','MarkerFaceAlpha',0.2);
        hold off
        
        [r1 p1]=corr(peaks_pfwe(ris_sig),abs(peakLags_pfwe(ris_sig)),'type','spearman','tail','left');
        [r2 p2]=corr(peaks(ris_unsig),abs(peakLags(ris_unsig)),'type','spearman','tail','left');
        
        legend({sprintf('R=%.2f (N=%d); p=%.3f',r1,length(ris_sig),p1),...
            sprintf('R=%.2f (N=%d); p=%.3f',r2,length(ris_unsig),p2)},'fontsize',12,'Location','northoutside');
        legend boxoff
        title({strrep([upper(exp(1)) exp(2:end)],'_',' ')});
        xlabel('Peak lag (sec)');
        ylabel('Peak R (z)');
        set(gca,'fontsize',14)
        xlim([-30 30])
        ylim([0 0.45])
        grid on
        
    end
end