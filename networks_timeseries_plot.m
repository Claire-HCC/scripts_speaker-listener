clear all
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
%
%  froidir='restFc_isc15PercMasked_50Overlap_cluster6';
%  networks={'AUD','vLAN','dLAN','DMNa','DMNb'};
% cols=jet(7);
% cols=cols([1 3 4 5 7],:);

froidir='restFc_isc30PercMasked_75Overlap_cluster6';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};
% cols=jet(7);
% cols=cols([1 3 4 5 6 7],:);
% cols(4,:)=[1 0.85 0];
cols=[74 37 170;
   66 152 181;
   40 71 52;
   219 11 91;
   221 135 141;
   115 83 29]./255;

crop_start=25;
crop_end=20;

for ei=[1 ];%[1 4 11 12 9 10];
    exp=exp_parameters.experiments{ei};
    
    figure('unit','centimeter','position',[0 0 18 6.5]);
    for ni=1:length(networks);
        %  network=networks{ni};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' networks{ni} ],'gdata','keptvox');
        
        [~,tn,listenerN]=size(gdata);
        keptT=(crop_start+1):(tn-crop_end);
        
        gdata(:,:,exp_parameters.subjects_excluded{ei})=NaN;
        gdata=nanmean(gdata,1);
        gdata=zscore(gdata(:,keptT,:),0,2);
        
        % ciplot_claire(squeeze(gdata)',keptT*exp_parameters.tr(ei),cols(ni,:),0.1);
        plot(keptT*exp_parameters.tr(ei),nanmean(squeeze(gdata)'),'color',cols(ni,:),'linewidth',4)
        hold on;
    end
end
% legend(networks,'orientation','horizontal','location','northoutside');
% legend boxoff;

% xlim([min(keptT) max(keptT)]);
xlim([320 520])
ylim([-1.2 1.2])
% set(gca,'xtick',[320:20:520])
% set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)

% line([0 0 ],get(gca,'ylim'),'color',gray);
% line(get(gca,'xlim'),[0 0],'color',gray);
% title([upper(exp(1)) strrep(exp(2:end),'_',' ')],'color','k');
xlabel('Time (sec)');
ylabel({'Group BOLD responses','to story (normalized)'});
set(gca,'fontsize',14);
hold off

rectangle('Position',[400 -1.2 20 2.4] ,'edgecolor','k','linewidth',3);
rectangle('Position',[500 -1.2 20 2.4] ,'edgecolor','k','linewidth',3);


xlim([400 420])
set(gcf,'position',[0 0 5 10]);
title('')
ylabel('')
xlabel('')
set(gca,'xtick',[],'ytick',[])

% xlim([500 520])