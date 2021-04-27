% close all
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=(roi_table.network);
rnames_complete=table2array(roi_table(:,3));
lags_tested={-20:20};
network_newOrder={
    'Auditory_Language',...
    'DMN2',...
    'Attention',...
    'Executive',...
    'DMN1',...
    'Visual'...
    };

binSize=30;
eis=[ 1 2 4 11 12 13 ];
fsize=[48 7]
figure('unit','centimeter','position',[0 0 fsize]);


for i=[1:length(eis)];
    ei=eis(i);
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/rois2rois_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaksMax' ],'rnames',...
            'rz','rzm','lags','keptT','p','pfwe','peaks_pfwe','peakLags_pfwe','pfdr','peaks_p','peakLags_p','peakLags','peaks');
        
        ris=[];
        for ni=1:length(network_newOrder);
            ris=[ris ; find(ismember(rnames,roi_table.region(ismember(roi_table.network,network_newOrder{ni}))))];
        end
        
        ris=ris(sum(~isnan(peakLags(ris,ris)),2)~=0);
          
        mask=zeros(size(peaks));
        mask(~isnan(peakLags_pfwe) )=1;
        
        subplot(1,6,i);
        
        imAlpha=ones(length(ris),length(ris));
        imAlpha(mask(ris,ris)==0)=0;
        
        peakLags_bined=discretize(peakLags,[-20 -15 -10 -6 -3 0 1 4  7  11 16 20]);
        im=imagesc(peakLags_bined(ris,ris),[2 10]);
        % im=imagesc(peaks(ris,ris),[-0.2 0.2]);
        
        %   im=imagesc(peakLags(ris,ris),[-10 10]);
        colormap jet
        set(im,'AlphaData',imAlpha)
        set(gca,'xtick',[],'ytick',[],'color',[0 0 0]);
        ax=gca;
        ax.XAxis.Visible='off';
        ax.XAxis.Label.Visible='on';
        ax.YAxis.Visible='off';
        ax.YAxis.Label.Visible='on';
        
        for ni=1:length(network_newOrder);
            network=network_newOrder{ni};
            temp=find(ismember(roi_table.network(ris),network));
            
            pos=[min(temp)-0.5 min(temp)-0.5 length(temp) length(temp)];
            rectangle('Position',pos ,'edgecolor','w','linewidth',0.5);
            
            if strcmp(network,'Auditory_Language'); network={'Auditory/','Language'};end
            text(pos(1)+pos(3)/2,pos(2)+pos(4)*0.5,network,'HorizontalAlignment','center','color','w','fontweight','bold','fontsize',9);
            
        end
        title(strrep([upper(exp(1)) exp(2:end)],'_',' '));
        ylabel('seed network');
        xlabel('target network');
        set(gca,'fontsize',14)
        
    end
end