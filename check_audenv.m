% clear all
% close all
% clear all
% set_parameters;
% timeUnit='tr' ;
% froidir='mor';
% load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
% rnames=table2array(roi_table(:,3));
% lags=-10:10;
% type='env';
% role='listener';
% rname='HG_L'
% ri=find(ismember(rnames,rname));
%
% fsize=[30,17];
% figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
%
% for ei=1:4;%3:4%:2;
%     exp=experiments{ei};
%
%     f=ls([expdir exp '/sound/' exp '_listener_cropped_audenv.mat' ]);
%     load([expdir experiments{ei} '/sound/' f ]);
%     load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
%     load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat' ],'data');
%     gdata(:,:,subjects_excluded{ei})=NaN;
%     gdata=nanmean(gdata,1);
%     data=nanmean(data,1);
%     [~,tn,listenerN]=size(gdata);
%     tn=min(min(tn,length(aud)),size(data,2));
%     r=[];
%     r=(lagcorr_claire(squeeze(gdata(:,:,:)),repmat(aud(1:tn),1,listenerN),lags))';
%     rs=(lagcorr(data(:,1:tn)',aud(1:tn),lags))';
%
%     subplot(2,4,ei);
%     plot(lags,r);
%     xlabel('lags (TR)');
%     ylabel('R');
%     title({rname,exp});
%     set(gca,'fontsize',14);
%     ylim([-0.4 0.4]);
%     grid on
%     xlim([-10 10])
%     set(gca,'xtick',-10:5:10);
%     grid on
%
%     subplot(2,4,ei+4);
%     plot(lags,rs,'r','linewidth',2);
%     hold on;
%     ciplot_claire(r,lags,'k',0.3);
%     hold off;
%     xlabel('lags (TR)');
%     ylabel('R');
%     title({rname,exp});
%     set(gca,'fontsize',14);
%     ylim([-0.4 0.4]);
%     grid on
%     xlim([-10 10])
%     set(gca,'xtick',-10:5:10)
% end
% legend('speaker','listener');

close all
clear all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
lags=-10:10;
type='env';
role='listener';
crop_start=0;
crop_end=0;
eis=[1 3 7 8];% 1:2;%1:4;

  %  figure;
    
for i=1:length(eis);
    ei=eis(i);

  %  subplot(2,4,i)
    exp=exp_parameters.experiments{ei};
    
    f=ls([expdir exp '/sound/*_listener_audenv.mat' ]);
    load([expdir exp '/sound/' f ]);
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_HG_L.mat' ],'gdata');
    gdata1=gdata;
    gdata1(:,:,exp_parameters.subjects_excluded{ei})=NaN;
    gdata1=nanmean(gdata1,1);
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_HG_R.mat' ],'gdata');
    gdata2=gdata;
    gdata2(:,:,exp_parameters.subjects_excluded{ei})=NaN;
    gdata2=nanmean(gdata2,1);
    
    [~,tn,listenerN]=size(gdata1);
    tn=min(min(tn,length(aud)));
    
    r1=(lagcorr_claire(aud(1:(tn-crop_start)),nanmean(gdata1((crop_start+1):tn),3)',lags))';
    r2=(lagcorr_claire(aud(1:(tn-crop_start)),nanmean(gdata2((crop_start+1):tn),3)',lags))';
    
%     plot(lags,r1,'linewidth',2);
%     hold on
%     plot(lags,r2,'linewidth',2);
%     hold off
%     legend('HG L','HG R')
%     xlabel('lags (TR)');
%     ylabel('R');
%     title(strrep(exp,'_hmhp',''));
%     set(gca,'fontsize',14);
%     
%     xlim([min(lags) max(lags)])
%     set(gca,'xtick',lags);
%     grid on
    
    gdata=(gdata1+gdata2)/2;
    r=(lagcorr_claire(repmat(aud(1:(tn-crop_start)),1,listenerN),squeeze(gdata(:,(crop_start+1):tn,:)),lags))';
    [peak lagi]=max(r,[],2);
    peakLags=lags(lagi);
    peakLags(exp_parameters.subjects_excluded{ei})=NaN;
    figure;
    subplot(1,2,1);
    hist(peakLags);
    title(exp);
    subplot(1,2,2);
     subjis=find(abs(peakLags-mode(peakLags))>5);
    if ~isempty(subjis);
        plot(lags,squeeze(r(subjis,:)));
      legend(num2str(subjis'));
    end
    %
    %     for si=1:listenerN;
    %
    %         plot(lags,r1(si,:),'linewidth',2);
    %                 hold on
    %         plot(lags,r2(si,:),'linewidth',2);
    %         hold off
    %         legend('HG L','HG R')
    %         xlabel('lags (TR)');
    %         ylabel('R');
    %         title({exp,num2str(si)});
    %         set(gca,'fontsize',14);
    %         %   ylim([-0.4 0.4]);
    %         xlim([-10 10])
    %         set(gca,'xtick',-10:5:10);
    %         grid on
    %         pause
    %     end
end


