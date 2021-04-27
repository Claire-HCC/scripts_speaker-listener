% close all
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=(roi_table.network);
rnames_complete=table2array(roi_table(:,3));
 roi_id=[(roi_table.id{:})]';
lags_tested={-15:15};
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


for i=1;%[1:length(eis)];
    ei=eis(i);
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs/LL_gg/voxs2voxs_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaksNearest' ],'keptvox_seed',...
            'lags','keptT','peakLags_pfwe','peakLags','peaks','peaks_pfwe');
        keptvox=keptvox_seed;
        
        vis=[];
        for ni=1:length(network_newOrder);
            network=networks{ni};
            fr = sprintf('%s/roi_mask/%s/mat/%s',expdir,froidir,network);
            load(fr,'roimask');
            vis=[vis find(roimask(keptvox)==1)'];
        end
        

        mask=zeros(size(peaks));
        mask(~isnan(peakLags_pfwe) )=1;
        
        subplot(1,6,i);
        
        imAlpha=ones(length(vis),length(vis));
        imAlpha(mask(vis,vis)==0)=0;
        
        peakLags_bined=discretize(peakLags,[-20 -15 -10 -6 -3 0 1 4  7  11 16 20]);
        im=imagesc(peakLags_bined(vis,vis),[2 10]);
        % im=imagesc(peaks(vis,vis),[-0.2 0.2]);
        
        %   im=imagesc(peakLags(vis,vis),[-10 10]);
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
            temp=find(ismember(roi_table.network(vis),network));
            
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