close all
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table','roi_newOrder');
networks=unique(roi_table.network);
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -40:40};
network_newOrder={
    'Auditory_Language',...
    'DMN2',...
    'DMN1',...
    'Attention',...
    'Executive',...
    'Visual'...
    };
roi_newOrder=[];
for ni=1:6;
    roi_newOrder=[roi_newOrder; find(ismember(roi_table.network,network_newOrder{ni}))];
end

for ei=3:4;%1:4;
    exp=experiments{ei};
    
    fsize=[7 16];
figure('unit','centimeter','position',[0 0 fsize]);

    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        fs= cellstr(ls([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/rois2rois*groupT*']));
        
        for fi=1:length(fs);
            load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/' fs{fi} ],...
                'rnames','t','r','rzm','lags','keptT','p','pfwe','peaks_pfwe','peakLags_pfwe','pfdr','peaks_pfdr','peakLags_pfdr','peakLags','peaks');
            
                subplot(2,1,fi)
            mask=zeros(size(peaks));
            mask(~isnan(peaks_pfdr) & peaks_pfdr>0)=1;
            %  mask(~isnan(pfdr(:,:,lags==0)) )=1;
            ris=roi_newOrder(ismember(roi_newOrder,find((nansum(mask,1))~=0)));
            
            imAlpha=ones(length(ris),length(ris));
            imAlpha(mask(ris,ris)==0)=0;
            
            im=imagesc(peakLags_pfdr(ris,ris),[-10 10]);
          
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
                
                pos=[min(temp)-0.5 min(temp)-0.5 length(temp)-0.1 length(temp)-0.1];
                rectangle('Position',pos ,'edgecolor',[0.3 0.3 0.3],'linewidth',1);
                
                if strcmp(network,'Auditory_Language'); network={'Auditory/','Language'};
                    text(pos(1)+pos(3)/2+2.5,pos(2)+pos(4)/2,network,'HorizontalAlignment','center','color',[0.3 0.3 0.3],'fontweight','bold','fontsize',9,'FontName', 'Arial');
                else
                    text(pos(1)+pos(3)/2,pos(2)+pos(4)/2,network,'HorizontalAlignment','center','color',[0.3 0.3 0.3],'fontweight','bold','fontsize',9,'FontName', 'Arial');
                end
            end
            title([upper(exp(1)) exp(2:end)]);
            %  title('Peak lag');
            ylabel('seed ROI');
            xlabel('target ROI');
            set(gca,'fontsize',14)
            
        end
    end
end