clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rgb=roi_table.rgb;
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
rgb=roi_table.rgb;
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
eis=[1 2 4 11 12 13];%
fsize=[45 7];
figure('unit','centimeter','position',[0 0 fsize],'color','k');

for i=1:length(eis);
    ei=eis(i);
    exp=experiments{ei};
    
    subplot(1,6,i);
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_selfself/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaksMax' ],...
        'networks','rzm','rz','lags','keptT','pfwe','peaks_pfwe','peakLags_pfwe','pfdr','peaks_p','peakLags_p','peakLags','peaks');
    
    for ni=1:6;
        seed=network_newOrder{ni};
        sdi_old=find(ismember(networks,seed));
        nis=cellfun(@(x) find(ismember(networks,x)),network_newOrder);
        
        line([min(lags) max(lags)],[0 0],'color',gray)
        hold on;
        
        temp=squeeze(rz(sdi_old,sdi_old,:,:))';
        
        ciplot_claire(temp,lags*1.5,[cols(ni+1,:) ],0.1);
        hold on;
        
        set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)
    %     ylim([-0.2 1])
        % xlim([min(lags) max(lags)]);
        xlim([-3 3])
          title(strrep(exp,'_',' '),'fontsize',14,'color','w');
     %   set(gca,'xtick',[-20:1:20],'xticklabel',[-20:1:20]*tr(ei));
        xlabel('Lag (sec)');
        ylabel('R (z)')
        set(gca,'fontsize',10)
        % grid on

    end
   
end
