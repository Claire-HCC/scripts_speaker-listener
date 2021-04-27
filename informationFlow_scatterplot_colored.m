%% cross stories LL bined consistency, networks2networks
clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='restFc_isc30PercMasked_50Overlap_cluster6';

lags=-20:20;

eis=1:5;
for i=1:length(eis);
    ei=eis(i);
    exp=exp_parameters.experiments{ei};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peaks' ],...
        'networks','rz','rzm','lags','keptT','pfdr_peaks','p_peaks','peakLags','peaks','pfdr_npeaks');
  
    peakLags_ll(:,ei)=peakLags(:);
    peaks_ll(:,ei)=peaks(:);
    sigs(:,ei)=(pfdr_peaks(:)<(.05) & ~isnan(p_peaks(:)) & (pfdr_npeaks(:)>.05 | isnan(pfdr_npeaks(:)) ));
    
end

fsize=[28 30];
figure('unit','centimeter','position',[0 0 fsize]);
subploti=reshape(1:(length(eis)*(length(eis)-1)),length(eis),length(eis)-1)';
for i1=1:length(eis);
    ei1=eis(i1);
    exp1=exp_parameters.experiments{ei1};
    
    for i2=1:length(eis);
        ei2=eis(i2);
        exp2=exp_parameters.experiments{ei2};
        
        if i2>i1;
            
            x=peakLags_ll(:,ei1);
            y=peakLags_ll(:,ei2);
                 
            shared=(sum(isnan(peakLags_ll(:,[ei1 ei2])),2)==0);
            sig=(sum(sigs(:,[ei1 ei2])==0,2)==0);
            
      %      triui=triu(ones(sqrt(length(x)),sqrt(length(x))),1);
       %     triui=triui(:);
            
           
           subplot(length(eis)-1,length(eis),subploti(i1,i2)) 
     %     subplot(4,5,i2-1) 
            hold on;
            
            t = 2*pi*rand(length(x),1);
            r = 0.5*sqrt(rand(length(x),1));
            x_jittered = x + r.*cos(t);
            y_jittered = y + r.*sin(t);
            
             if sum(sig)~=0;
              scatter(x_jittered(sig),y_jittered(sig),20,'k','filled','MarkerFaceAlpha',1);
             end
             if sum(sig)>1;
                   [r1 p1]=corr(x(sig),y(sig ),'type','spearman','tail','right');
             end
            
            [r2 p2]=corr(x(shared),y(shared ),'type','spearman','tail','right');
         
             scatter(x_jittered(shared),y_jittered(shared),20,'k','filled','MarkerFaceAlpha',0.2);

            set(gca,'xtick',[-20:4:20],'xticklabel',[-20:4:20]*exp_parameters.tr(ei));
            set(gca,'ytick',[-20:4:20],'xticklabel',[-20:4:20]*exp_parameters.tr(ei));
            xlim([min(x)*0.9 max(x)*1.1]);
            ylim([min(y)*0.9 max(y)*1.1]);
            
            yticklabels( get(gca,'ytick')*exp_parameters.tr(ei2));
            set(gca,'fontsize',12);
         
            legend({sprintf('R=%.2f (N=%d)',r1,sum(sig )),...
                sprintf('R=%.2f (N=%d)',r2,sum(shared ))},'fontsize',8,'Location','northoutside');
            legend boxoff
            xlabel({[upper(exp1(1)) strrep(exp1(2:end),'_',' ') ]},'fontsize',12,'fontweight','bold');
            ylabel({[upper(exp2(1)) strrep(exp2(2:end),'_',' ') ]},'fontsize',12,'fontweight','bold');
        end
    end
end

%
% %% cross stories LL consistency
% clear all
% close all
% set_parameters;
% timeUnit='tr' ;
% froidir='mor';
% load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
% rgb=roi_table.rgb;
% rnames=table2array(roi_table(:,3));
% roi_ids=cell2mat(roi_table.id);
% rgb=roi_table.rgb;
% lags=-10:10;
%
% seed='vPCUN'
% for ei=1:4;
%     exp=experiments{ei};
%     load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT_peaks' ],...
%         'rnames','t','r','rzg','lags','keptT','p','pfdr','peaks_pfdr','peakLags_pfdr','pfdr','peaks_pfdr','peakLags_pfdr','peakLags','peaks');
%
%     peakLags_ll_pfdr(:,ei)=peakLags_pfdr;
%     peaks_ll_pfdr(:,ei)=peaks_pfdr;
%     peakLags_ll(:,ei)=peakLags;
%     peaks_ll(:,ei)=peaks;
% end
% peakLags_ll_pfdr(peaks_ll_pfdr<0)=NaN;
% peakLags_ll(peaks_ll<0)=NaN;
%
% % to reverse the direction of information flow in figures.
% % peakLags_ll_pfdr=-peakLags_ll_pfdr;
% % peakLags_ll=-peakLags_ll;
%
% fsize=[10.5 11];
% subploti=reshape([1:12],4,3);
% for ei1=4;%:4;
%     exp1=experiments{ei1};
%
%     for ei2=3;%1:4;
%         exp2=experiments{ei2};
%
%         if ei1>ei2;
%             ris=find(sum(isnan(peakLags_ll_pfdr(:,[ei1 ei2])),2)==0);
%             ris_unsig=find(~ismember(1:length(rnames),ris)' & (sum(isnan(peakLags_ll(:,[ei1 ei2])),2)==0));
%             n=length(ris);
%             x=peakLags_ll(:,ei1);
%             y=peakLags_ll(:,ei2);
%             [r1 p1]=corr(x(ris),y(ris),'type','spearman','tail','right');
%             [r2 p2]=corr(x([ris; ris_unsig]),y([ris; ris_unsig]),'type','spearman','tail','right');
%
%
%             %   figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize]);
%             % subplot(3,4,subploti(ei1,ei2));
%             fsize=[10.5 11];
%             figure('unit','centimeter','position',[0 0 fsize]);
%             x_jittered=x+(rand(length(x),1)-0.5)*0.5;
%             y_jittered=y+(rand(length(y),1)-0.5)*0.5;
%
%             scatter(x_jittered(ris_unsig),y_jittered(ris_unsig),50,'k','filled','MarkerFaceAlpha',.3);
%             hold on
%             % scatter(x_jittered(ris),y_jittered(ris),16,'r','filled','MarkerFaceAlpha',.2,'markerEdgeColor','r');
%             scatter(x_jittered(ris),y_jittered(ris),80,rgb(ris,:),'linewidth',1.5);
%
%             xticklabels( get(gca,'xtick')*tr(ei));
%             yticklabels( get(gca,'ytick')*tr(ei));
%
%             xlim([min(x)-1 max(x)+1]);
%             ylim([min(y)-1 max(y)+1]);
%             set(gca,'fontsize',15)
%             title({[strrep(seed,'_',' ') ' seed'],sprintf('Spearman''s R=%.2f (N=%d overlapped ROIs)',r1,length(ris)),...
%                 sprintf('Spearman''s R=%.2f (N=%d all ROIs)',r2,length([ris; ris_unsig]))},'fontsize',12);
%             xlabel({[upper(exp1(1)) exp1(2:end) ' Peak Lag (sec)']},'fontsize',16,'fontweight','bold');
%             ylabel({[upper(exp2(1)) exp2(2:end) ' Peak Lag (sec)']},'fontsize',16,'fontweight','bold');
%         end
%     end
% end
%
%
%
% %% cross seeds
% clear all
%
% set_parameters;
% timeUnit='tr' ;
% froidir='mor';
% load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
% rnames=table2array(roi_table(:,3));
% roi_ids=cell2mat(roi_table.id);
% rgb=roi_table.rgb;
% lags=-10:10;
% seeds={'vPCUN','pANG_L','HG_L'};
% for ei=3:4;
%     exp=experiments{ei};
%
%     for sdi=1:length(seeds);
%         seed=seeds{sdi};
%
%         load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT_peaks' ],...
%             'rnames','t','r','rzg','lags','keptT','p','pfdr','peaks_pfdr','peakLags_pfdr','pfdr','peaks_pfdr','peakLags_pfdr','peakLags','peaks');
%
%         peakLags_ll(:,ei,sdi)=peakLags;
%         peaks_ll(:,ei,sdi)=peaks;
%         peakLags_ll_pfdr(:,ei,sdi)=peakLags_pfdr;
%         peaks_ll_pfdr(:,ei,sdi)=peaks_pfdr;
%     end
% end
% peakLags_ll(peaks_ll<0)=NaN;
% peakLags_ll_pfdr(peaks_ll_pfdr<0)=NaN;
% % to reverse the direction of information flow in figures.
% % peakLags_ll=-peakLags_ll;
% % peakLags_ll_pfdr=-peakLags_ll_pfdr;
%
% for ei=3:4;
%     exp=experiments{ei};
%
%     for sdi1=1:length(seeds);
%         for sdi2=1:length(seeds);
%             if sdi1>sdi2;
%
%                 ris=find(sum(isnan(peakLags_ll_pfdr(:,ei,[sdi1 sdi2])),3)==0);
%                 ris_unsig=find(~ismember(1:length(rnames),ris)' & sum(isnan(peakLags_ll(:,ei,[sdi1 sdi2])),3)==0);
%
%                 n=length(ris);
%                 x=squeeze(peakLags_ll(:,ei,sdi1));
%                 y=squeeze(peakLags_ll(:,ei,sdi2));
%                 [r1 p1]=corr(x(ris),y(ris),'type','spearman','tail','right');
%                 [r2 p2]=corr(x([ris; ris_unsig]),y([ris; ris_unsig]),'type','spearman','tail','right');
%
%                 fsize=[10.5 11];
%                 figure('unit','centimeter','position',[0 0 fsize]);
%                 x_jittered=x+(rand(length(x),1)-0.5)*0.8;
%                 y_jittered=y+(rand(length(y),1)-0.5)*0.8;
%                 scatter(x_jittered(ris_unsig),y_jittered(ris_unsig),50,'k','filled','MarkerFaceAlpha',.3);
%                 hold on
%                 scatter(x_jittered(ris),y_jittered(ris),60,rgb(ris,:),'linewidth',2);
%                 % scatter(x_jittered(ris),y_jittered(ris),16,'r','filled','MarkerFaceAlpha',.2,'markerEdgeColor','r');
%                 xlabel({[strrep(seeds{sdi1},'_',' ')],'before seed--Lag(sec)--after seed'});
%                 ylabel({'before seed--Lag(sec)--after seed',[strrep(seeds{sdi2},'_','')]});
%
%                 xlim([min(x)-1 max(x)+1]);
%                 ylim([min(y)-1 max(y)+1]);
%                 xticklabels( get(gca,'xtick')*tr(ei));
%                 yticklabels( get(gca,'ytick')*tr(ei));
%                 set(gca,'fontsize',16)
%                 title({[upper(exp(1)) exp(2:end)],sprintf('Spearman''s R=%.2f (N=%d overlapped ROIs)',r1,length(ris)),...
%                     sprintf('Spearman''s R=%.2f (N=%d all ROIs)',r2,length([ris; ris_unsig]))},'fontsize',12);
%             end
%         end
%     end
% end
%
% %% SL vs. LL
% clear all
%
% set_parameters;
% timeUnit='tr' ;
% froidir='mor';
% load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
% rnames=table2array(roi_table(:,3));
% roi_ids=cell2mat(roi_table.id);
% rgb=roi_table.rgb;
% lags=-10:10;
% seed='vPCUN';
% for ei=1:4;
%     exp=experiments{ei};
%
%     load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT_peaks' ],...
%         'rnames','t','r','rzg','lags','keptT','p','pfdr','peaks_pfdr','peakLags_pfdr','pfdr','peaks_pfdr','peakLags_pfdr','peaks','peakLags');
%
%     peakLags_ll(:,ei)=peakLags;
%     peaks_ll(:,ei)=peaks;
%
%     peakLags_ll_pfdr(:,ei)=peakLags_pfdr;
%     peaks_ll_pfdr(:,ei)=peaks_pfdr;
%
%     load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/SL_each/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT_peaks' ],...
%         'rnames','t','r','rzg','lags','keptT','p','pfdr','peaks_pfdr','peakLags_pfdr','pfdr','peaks_pfdr','peakLags_pfdr','peaks','peakLags');
%
%     peakLags_sl(:,ei)=peakLags;
%     peaks_sl(:,ei)=peaks;
%
%     peakLags_sl_pfdr(:,ei)=peakLags_pfdr;
%     peaks_sl_pfdr(:,ei)=peaks_pfdr;
% end
%
% % peakLags_sl=-peakLags_sl;
% % peakLags_sl_pfdr=-peakLags_sl_pfdr;
% % to reverse the direction of information flow in figures.
% % peakLags_ll=-peakLags_ll;
% % peakLags_ll_pfdr=-peakLags_ll_pfdr;
%
% for ei=3:4;
%     exp=experiments{ei};
%     ris=find(sum(isnan(peakLags_sl_pfdr(:,ei)),3)==0 & sum(isnan(peakLags_ll_pfdr(:,ei)),3)==0);
%     ris_unsig=find(~ismember(1:length(rnames),ris)' & sum(isnan(peakLags_sl(:,ei)),3)==0);
%
%     if length(ris)>=2;
%         n=length(ris);
%         x=squeeze(peakLags_ll(:,ei));
%         y=squeeze(peakLags_sl(:,ei));
%         [r1 p1]=corr(x(ris),y(ris),'type','spearman','tail','right');
%         [r2 p2]=corr(x([ris; ris_unsig]),y([ris; ris_unsig]),'type','spearman','tail','right');
%
%
%         fsize=[11 11.5];
%         figure('unit','centimeter','position',[0 0 fsize]);
%         x_jittered=x+(rand(length(x),1)-0.5)*0.8;
%         y_jittered=y+(rand(length(y),1)-0.5)*0.8;
%         scatter(x_jittered(ris_unsig),y_jittered(ris_unsig),50,'k','filled','MarkerFaceAlpha',.3);
%         hold on
%         scatter(x_jittered(ris),y_jittered(ris),60,rgb(ris,:),'linewidth',2);
%         % scatter(x_jittered(ris),y_jittered(ris),16,'r','filled','MarkerFaceAlpha',.2,'markerEdgeColor','r');
%         xlim([min(min(y),min(x))-1 max(max(y),max(x))+1]);
%         ylim([min(min(y),min(x))-1 max(max(y),max(x))+1]);
%         xticklabels( get(gca,'xtick')*tr(ei));
%         yticklabels( get(gca,'ytick')*tr(ei));
%         xlabel({['Listener-Listener Peak Lag'],'before seed--Lag(sec)--after seed'});
%         ylabel({'before seed--Lag(sec)--after seed','Speaker-Listener Peak Lag'});
%
%         set(gca,'fontsize',16)
%
%         title({[upper(exp(1)) exp(2:end)],sprintf('Spearman''s R=%.2f (N=%d overlapped ROIs)',r1,length(ris)),...
%             sprintf('Spearman''s R=%.2f (N=%d all ROIs)',r2,length([ris; ris_unsig]))},'fontsize',12);
%     end
% end
%
% %% SL  across stories
% clear all
% set_parameters;
% timeUnit='tr' ;
% froidir='mor';
% load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
% rnames=table2array(roi_table(:,3));
% roi_ids=cell2mat(roi_table.id);
% rgb=roi_table.rgb;
% lags=-10:10;
% seed='vPCUN'
% for ei=3:4;
%     exp=experiments{ei};
%
%     load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/SL_each/' seed '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT_peaks'],...
%         'rnames','t','r','rzg','lags','keptT','p','pfdr','peaks_pfdr','peakLags_pfdr','pfdr','peaks_pfdr','peakLags_pfdr','peakLags','peaks');
%
%     peakLags_sl_pfdr(:,ei)=peakLags_pfdr;
%     peaks_sl_pfdr(:,ei)=peaks_pfdr;
%
%     peakLags_sl(:,ei)=peakLags;
%     peaks_sl(:,ei)=peaks;
% end
%
% peakLags_sl(peaks_sl<0)=NaN;
% peakLags_sl_pfdr(peaks_sl_pfdr<0)=NaN;
%
% % to reverse the direction of information flow in figures.
% % peakLags_sl=-peakLags_sl;
% % peakLags_sl_pfdr=-peakLags_sl_pfdr;
%
% for ei1=4;%1:4;
%     exp1=experiments{ei1};
%
%     for ei2=3;%1:4;
%         exp2=experiments{ei2};
%
%         if ei1>ei2;
%             ris=find(sum(isnan(peakLags_sl_pfdr(:,[ei1 ei2])),2)==0);
%             ris_unsig=find(~ismember(1:length(rnames),ris)' & (sum(isnan(peakLags_sl(:,[ei1 ei2])),2)==0));
%
%             x=peakLags_sl(:,ei1);
%             y=peakLags_sl(:,ei2);
%             [r1 p1]=corr(x(ris),y(ris),'type','spearman','tail','right');
%             [r2 p2]=corr(x([ris; ris_unsig]),y([ris; ris_unsig]),'type','spearman','tail','right');
%
%             %   figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize]);
%             % subplot(3,4,subploti(ei1,ei2));
%             fsize=[10.5 11];
%             figure('unit','centimeter','position',[0 0 fsize]);
%             x_jittered=x+(rand(length(x),1)-0.5)*0.5;
%             y_jittered=y+(rand(length(y),1)-0.5)*0.5;
%
%             scatter(x_jittered(ris_unsig),y_jittered(ris_unsig),50,'k','filled','MarkerFaceAlpha',.3);
%             hold on
%             % scatter(x_jittered(ris),y_jittered(ris),16,'r','filled','MarkerFaceAlpha',.2,'markerEdgeColor','r');
%             scatter(x_jittered(ris),y_jittered(ris),80,rgb(ris,:),'linewidth',1.5);
%
%             xlim([min(x)-1 max(x)+1]);
%             ylim([min(y)-1 max(y)+1]);
%             xticklabels( get(gca,'xtick')*tr(ei));
%             yticklabels( get(gca,'ytick')*tr(ei));
%             xlabel({[upper(exp1(1)) exp1(2:end)],'before seed--Lag(sec)--after seed'});
%             ylabel({'before seed--Lag(sec)--after seed',[upper(exp2(1)) exp2(2:end)]});
%
%             set(gca,'fontsize',15)
%
%             title({[strrep(seed,'_',' ') ' seed'],sprintf('Spearman''s R=%.2f (N=%d overlapped ROIs)',r1,length(ris)),...
%                 sprintf('Spearman''s R=%.2f (N=%d all ROIs)',r2,length([ris; ris_unsig]))},'fontsize',12);
%             xlabel({[upper(exp1(1)) exp1(2:end) ' Peak Lag (sec)']},'fontsize',16,'fontweight','bold');
%             ylabel({[upper(exp2(1)) exp2(2:end) ' Peak Lag (sec)']},'fontsize',16,'fontweight','bold');
%         end
%     end
% end
%
% %% cross stories LL consistency, rois2rois
% clear all
% close all
% set_parameters;
% timeUnit='tr' ;
% froidir='mor';
% load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
% rgb=roi_table.rgb;
% rnames=table2array(roi_table(:,3));
% roi_ids=cell2mat(roi_table.id);
% rgb=roi_table.rgb;
% lags=-10:10;
%
% for ei=1:4;
%     exp=experiments{ei};
%     load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/rois2rois_lag' num2str(min(lags)) '-' num2str(max(lags)) '_groupT2tail_peaks' ],...
%         'rnames','t','r','rzm','lags','keptT','p','pfdr','peaks_pfdr','peakLags_pfdr','pfdr','peaks_pfdr','peakLags_pfdr','peaks_pfwe','peakLags_pfwe','peakLags','peaks');
%
%     peakLags_ll_pfdr(:,ei)=peakLags_pfdr(:);
%     peaks_ll_pfdr(:,ei)=peaks_pfdr(:);
%     peakLags_ll(:,ei)=peakLags(:);
%     peaks_ll(:,ei)=peaks(:);
%     peakLags_ll_pfwe(:,ei)=peakLags_pfwe(:);
%     peaks_ll_pfwe(:,ei)=peaks_pfwe(:);
%
% end
%
% % remove negative peaks
% peakLags_ll_pfwe(peaks_ll<0)=NaN;
% peakLags_ll_pfdr(peaks_ll<0)=NaN;
% peakLags_ll(peaks_ll<0)=NaN;
%
% fsize=[10.5 11];
% subploti=reshape([1:12],4,3);
% for ei1=4;
%     exp1=experiments{ei1};
%
%     for ei2=3;
%         exp2=experiments{ei2};
%
%         if ei1>ei2;
%
%             x=peakLags_ll_pfdr(:,ei1);
%             y=peakLags_ll_pfdr(:,ei2);
%             x_unsig=peakLags_ll(:,ei1);
%             y_unsig=peakLags_ll(:,ei2);
%             ris=find(sum(isnan([x y]),2)==0);
%             ris_unsig=find(~ismember(1:size(peakLags_ll,1),ris)' & (sum(isnan([x_unsig y_unsig]),2)==0));
%             n=length(ris);
%
%             triui=triu(ones(sqrt(length(x)),sqrt(length(x))),1);
%             triui=find(triui(:));
%
%             [r1 p1]=corr(x(intersect(ris , triui)),y(intersect(ris , triui)),'type','spearman','tail','right');
%             [r2 p2]=corr(x_unsig(intersect([ris; ris_unsig],triui(:))),y_unsig(intersect([ris; ris_unsig],triui(:))),'type','spearman','tail','right');
%
%
%             %   figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize]);
%             % subplot(3,4,subploti(ei1,ei2));
%             fsize=[10.5 11];
%             figure('unit','centimeter','position',[0 0 fsize]);
%
%             t = 2*pi*rand(length(x),1);
%             r = 0.5*sqrt(rand(length(x),1));
%             x_jittered = x + r.*cos(t);
%             y_jittered = y + r.*sin(t);
%             x_unsig_jittered = x_unsig + r.*cos(t);
%             y_unsig_jittered = y_unsig + r.*sin(t);
%
%             scatter(x_unsig_jittered(ris_unsig),y_unsig_jittered(ris_unsig),5,'k','filled','MarkerFaceAlpha',.2);
%             hold on
%
%             scatter(x_jittered(ris),y_jittered(ris),5,'k','filled','MarkerFaceAlpha',1);
%
%             xticklabels( get(gca,'xtick')*tr(ei));
%             yticklabels( get(gca,'ytick')*tr(ei));
%
%             xlim([min(x)-1 max(x)+1]);
%             ylim([min(y)-1 max(y)+1]);
%             set(gca,'fontsize',15)
%             title({sprintf('Spearman''s R=%.2f (N=%d overlapped ROI pairs)',r1,length(ris)),...
%                 sprintf('Spearman''s R=%.2f (N=%d all ROI pairs)',r2,length([ris; ris_unsig]))},'fontsize',10);
%             xlabel({[upper(exp1(1)) exp1(2:end) ' Peak Lag (sec)']},'fontsize',16,'fontweight','bold');
%             ylabel({[upper(exp2(1)) exp2(2:end) ' Peak Lag (sec)']},'fontsize',16,'fontweight','bold');
%         end
%     end
% end
%
% %% cross stories LL bined consistency, rois2rois, bined
% clear all
% close all
% set_parameters;
% timeUnit='tr' ;
% froidir='mor';
% load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
% rgb=roi_table.rgb;
% rnames=table2array(roi_table(:,3));
% roi_ids=cell2mat(roi_table.id);
% rgb=roi_table.rgb;
% lags=-10:10;
%
% for ei=3:4;
%     exp=experiments{ei};
%
%     fs= cellstr(ls([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/rois2rois*groupT*']));
%
%     for fi=1:length(fs);
%         load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out_bined/' fs{fi} ],...
%             'rnames','t','r','rzm','lags','keptT','p','pfwe','peaks_pfwe','peakLags_pfwe','pfdr','peaks_pfdr','peakLags_pfdr','peakLags','peaks');
%
%         peakLags_ll_pfdr(:,ei,fi)=peakLags_pfdr(:);
%         peaks_ll_pfdr(:,ei,fi)=peaks_pfdr(:);
%         peakLags_ll(:,ei,fi)=peakLags(:);
%         peaks_ll(:,ei,fi)=peaks(:);
%         peakLags_ll_pfwe(:,ei,fi)=peakLags_pfwe(:);
%         peaks_ll_pfwe(:,ei,fi)=peaks_pfwe(:);
%     end
% end
%
% % remove negative peaks
% peakLags_ll_pfwe(peaks_ll<0)=NaN;
% peakLags_ll_pfdr(peaks_ll<0)=NaN;
% peakLags_ll(peaks_ll<0)=NaN;
%
% binNames={'first','last'};
% fsize=[10.5 11];
% subploti=reshape([1:12],4,3);
% for ei1=3:4;
%     exp1=experiments{ei1};
%
%     for ei2=3:4;
%         exp2=experiments{ei2};
%
%         if ei1>=ei2;
%             for fi1=1:2;
%                 for fi2=1:2;
%                     if fi1>=fi2;
%
%                         x=peakLags_ll_pfdr(:,ei1,fi1);
%                         y=peakLags_ll_pfdr(:,ei2,fi2);
%                         x_unsig=peakLags_ll(:,ei1,fi1);
%                         y_unsig=peakLags_ll(:,ei2,fi2);
%                         ris=find(sum(isnan([x y]),2)==0);
%                         ris_unsig=find(~ismember(1:size(peakLags_ll,1),ris)' & (sum(isnan([x_unsig y_unsig]),2)==0));
%                         n=length(ris);
%
%                         triui=triu(ones(sqrt(length(x)),sqrt(length(x))),1);
%                         triui=find(triui(:));
%
%                         [r1 p1]=corr(x(intersect(ris , triui)),y(intersect(ris , triui)),'type','spearman','tail','right');
%                         [r2 p2]=corr(x_unsig(intersect([ris; ris_unsig],triui(:))),y_unsig(intersect([ris; ris_unsig],triui(:))),'type','spearman','tail','right');
%
%                         %   figure('unit','centimeter','position',[0 0 fsize],'papersize',fsize,'paperposition',[0 0 fsize]);
%                         % subplot(3,4,subploti(ei1,ei2));
%                         fsize=[10.5 11];
%                         figure('unit','centimeter','position',[0 0 fsize]);
%
%                         t = 2*pi*rand(length(x),1);
%                         r = 0.5*sqrt(rand(length(x),1));
%                         x_jittered = x + r.*cos(t);
%                         y_jittered = y + r.*sin(t);
%                         x_unsig_jittered = x_unsig + r.*cos(t);
%                         y_unsig_jittered = y_unsig + r.*sin(t);
%
%                         scatter(x_unsig_jittered(ris_unsig),y_unsig_jittered(ris_unsig),5,'k','filled','MarkerFaceAlpha',.2);
%                         hold on
%
%                         scatter(x_jittered(ris),y_jittered(ris),5,'k','filled','MarkerFaceAlpha',1);
%
%                         xticklabels( get(gca,'xtick')*tr(ei));
%                         yticklabels( get(gca,'ytick')*tr(ei));
%
%                         xlim([min(x)-1 max(x)+1]);
%                         ylim([min(y)-1 max(y)+1]);
%                         set(gca,'fontsize',15)
%                         title({sprintf('Spearman''s R=%.2f (N=%d overlapped ROI pairs)',r1,length(ris)),...
%                             sprintf('Spearman''s R=%.2f (N=%d all ROI pairs)',r2,length([ris; ris_unsig]))},'fontsize',10);
%                         xlabel({[upper(exp1(1)) exp1(2:end) ' Peak Lag (sec)'],[binNames{fi1} ' 250 TR']},'fontsize',16,'fontweight','bold');
%                         ylabel({[upper(exp2(1)) exp2(2:end) ' Peak Lag (sec)'],[binNames{fi2} ' 250 TR']},'fontsize',16,'fontweight','bold');
%                     end
%                 end
%             end
%         end
%     end
% end
%
%
%
% %% cross stories LL bined consistency, networks2networks, pca
% clear all
% close all
% set_parameters;
% timeUnit='tr' ;
% froidir='mor';
% load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
% rgb=roi_table.rgb;
% rnames=table2array(roi_table(:,3));
% roi_ids=cell2mat(roi_table.id);
% rgb=roi_table.rgb;
% lags=-10:10;
% binSize=30;
%
% eis=[1:2 4 11 12 9 10 ];% 1:2;%1:4;
% for i=1:length(eis);
%     ei=eis(i);
%     exp=experiments{ei};
%     pci=1;
%     %  if ei==1; pci=2; else; pci=1; end
%     load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out_bined/networks2networks' '_binSize' num2str(binSize)  '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_pc' num2str(pci)],...
%         'pc','peaks','peakLags');
%     peakLags_ll(:,ei)=peakLags(:);
%     peaks_ll(:,ei)=peaks(:);
% end
%
% fsize=[30 28];
% figure('unit','centimeter','position',[0 0 fsize]);
% subploti=reshape(1:(length(eis)*(length(eis)-1)),length(eis),length(eis)-1)';
% for i1=1:length(eis);
%     ei1=eis(i1);
%     exp1=experiments{ei1};
%
%     for i2=1:length(eis);
%         ei2=eis(i2);
%         exp2=experiments{ei2};
%
%         if i2>i1;
%
%             x=peakLags_ll(:,ei1);
%             y=peakLags_ll(:,ei2);
%
%             ris=find(sum(isnan([x y]),2)==0);
%
%             n=length(ris);
%
%             triui=triu(ones(sqrt(length(x)),sqrt(length(x))),1);
%             triui=find(triui(:));
%
%             [r1 p1]=corr(x(intersect(ris , triui)),y(intersect(ris , triui)),'type','spearman','tail','right');
%             subplot(length(eis)-1,length(eis),subploti(i1,i2))
%
%             t = 2*pi*rand(length(x),1);
%             r = 0.5*sqrt(rand(length(x),1));
%             x_jittered = x + r.*cos(t);
%             y_jittered = y + r.*sin(t);
%
%             scatter(x_jittered(ris),y_jittered(ris),15,'k','filled','MarkerFaceAlpha',1);
%
%             xticklabels( get(gca,'xtick')*tr(ei1));
%             yticklabels( get(gca,'ytick')*tr(ei2));
%
%             xlim([min(x)-1 max(x)+1]);
%             ylim([min(y)-1 max(y)+1]);
%             set(gca,'fontsize',12)
%             title(sprintf('R=%.2f, p=%.2f',r1,p1),'fontsize',12);
%             xlabel({[upper(exp1(1)) strrep(exp1(2:end),'_',' ') ]},'fontsize',10,'fontweight','bold');
%             ylabel({[upper(exp2(1)) strrep(exp2(2:end),'_',' ') ]},'fontsize',10,'fontweight','bold');
%         end
%     end
% end