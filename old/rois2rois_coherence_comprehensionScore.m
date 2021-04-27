clear all
% close all
set_parameters;
win=66; %(100/1.5);
eis=[1 2 4 11 12 13];
timeUnit='tr';
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
names=table2array(roi_table(:,3));
crop_start=25;
crop_end=20;

Fs=1/1.5;
rnames_selected={'vPCUN','HG_L'};
cols=[0 0 1;1 0 0];

eis=[1 2 11 12 9 10];
fsize=[30 18];
figure('unit','centimeter','position',[0 0 fsize]);
cols=[0 0 1; 1 0 0];
% fbs=
% for eii=1:length(eis);
%     ei=eis(eii);
%     exp=experiments{ei};
%     load([expdir '/' exp '/bhv/comprehensionScore.mat'  ],'score');
%     load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/rois2rois_coherence'  ],'pha','coh','keptT','freq','Cxy','rnames_selected');
%     subplot(2,3,eii);
%     
%     subjs=find(~isnan(squeeze(coh(1,1,1,:))));
%     for ri=1:length(rnames_selected);
%         scatter(squeeze(sum(coh(ri,ri,freq>0 & freq<0.1,subjs),3)),score(subjs),40,cols(ri,:),'filled');
%         [r(ei,ri,1), p(ei,ri,1)]= corr(squeeze(sum(coh(ri,ri,freq>0 & freq<0.1,subjs),3)),score(subjs));
%         hold on
%         scatter(squeeze(sum(coh(ri,ri,freq>0.1 & freq<0.2,subjs),3)),score(subjs),40,cols(ri,:))
%         [r(ei,ri,2), p(ei,ri,2,2)]= corr(squeeze(sum(coh(ri,ri,freq>0.1 & freq<0.2,subjs),3)),score(subjs));
%         
%     end
%     
%     title([upper(exp(1)) strrep(exp(2:end),'_',' ')],'color','k');
%     
%     hold off
% end
eis=[1 2 11 12 9 10];
cols=zeros(61,3);
cols(1,:)=[1 0 0];
cols(11,:)=[0 0 1];
figure
for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    load([expdir '/' exp '/bhv/comprehensionScore.mat'  ],'score');
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/rois2rois_coherence'  ],'pha','coh','keptT','freq','Cxy','rnames_selected');
    subplot(2,4,eii);
    
    subjs=find(~isnan(squeeze(coh(1,1,1,:))));
    for ri=[1 11]
        for fi=1:length(freq);
            [r(ri,fi,eii) p(ri,fi,eii)]= corr(squeeze(coh(ri,ri,fi,subjs)),score(subjs),'type','spearman');
        end
        pfwe(ri,:,eii)=p(ri,:,eii)*length(freq);
        [~,~,pfdr(ri,:,eii)]=fdr(p(ri,:,eii));
        
        pc=plot(freq(2:end),squeeze(r(ri,2:end,eii)),'linewidth',2,'color',cols(ri,:));
        hold on
        sig=[squeeze(pfdr(ri,:,eii))<.05];
        
        if sum(sig)>0;
        scatter(freq(sig),r(ri,sig,eii),40,'markerfacecolor',get(pc,'color'));
        end
    end
    ylim([-0.6 0.6]);
    xlim([min(freq) max(freq)]);
end

figure;
ciplot_claire(atanh(squeeze(r(11,:,:)))',freq,'b',0.3)
hold on;
ciplot_claire(atanh(squeeze(r(1,:,:)))',freq,'r',0.3)
grid on;
xlabel('Frequency (Hz)')
set(gca,'fontsize',14)
ylabel({'Correlation between','coherence and comprehensino score'})


