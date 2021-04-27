clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
crop_start=25;
crop_end=20;
cutoff=[0.1];
Fs=1/1.5;
rnames_selected={'HG_L','vPCUN'}

for ei=[1 ]
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_iscmasked_HG_L.mat' ],'gdata');
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
    
    gdata_lp=[];
    gdata_hp=[];
    gdata_orig=[];
    for ri=1:length(rnames_selected);
        rname=rnames_selected{ri};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_iscmasked_'  rname '.mat' ],'gdata');

        gdata(:,:,subjects_excluded{ei})=NaN;
        gdata=nanmean(gdata,1);
        gdata=zscore(gdata(:,keptT,:),0,2);
        gdata_orig(ri,:,:)=gdata;
        for si=1:listenerN;
            gdata_lp(ri,:,si)=filter1('lp',double(gdata(:,:,si)),'fc',cutoff,'fs',Fs);
            gdata_hp(ri,:,si)=filter1('hp',double(gdata(:,:,si)),'fc',cutoff,'fs',Fs);
        end
    end
end



fsize=[30 22];
figure('unit','centimeter','position',[0 0 fsize],'color','k');
cols=jet(7);
cols=cols([2 6],:);
gray=[0.7 0.7 0.7];
for tgi=[1 2];
    
    subplot(3,2,2);
    plot(lags,circularlagcorr(squeeze(nanmean(gdata_orig(1,:,1:9),3))',squeeze(nanmean(gdata_orig(tgi,:,10:18),3))',lags),'linewidth',2);
    hold on
    set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray);
    xlabel('Lag (TR=1.5 sec)');
    ylabel('R')
    ylim([-0.5 0.5]);
    title([rnames_selected{1} ' seed'])
    
    subplot(3,2,1);
    ciplot_claire(squeeze(gdata_orig(tgi,:,:))',keptT,cols(tgi,:),0.3);
    hold on
end
xlim([min(keptT) max(keptT)]);
ylim([-1.3 1.3])

set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)

line([0 0 ],get(gca,'ylim'),'color',gray);
line(get(gca,'xlim'),[0 0],'color',gray);
title([upper(exp(1)) strrep(exp(2:end),'_',' ') ', original'],'color','w');
xlabel('Time (TR)');
ylabel(['fMRI signal']);
set(gca,'fontsize',14);
grid on
hold off


for tgi=[1 2];
    subplot(3,2,4);
    plot(lags,circularlagcorr(squeeze(nanmean(gdata_hp(1,:,1:9),3))',squeeze(nanmean(gdata_hp(tgi,:,10:18),3))',lags),'linewidth',2);
    hold on
    set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)
    ylim([-0.5 0.5]);
    
    subplot(3,2,3);
    ciplot_claire(squeeze(gdata_hp(tgi,:,:))',keptT,cols(tgi,:),0.3);
    
    hold on
end
%xlim([50 350]);
ylim([-1 1])
xlim([min(keptT) max(keptT)]);
gray=[0.7 0.7 0.7];
set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)

line([0 0 ],get(gca,'ylim'),'color',gray);
line(get(gca,'xlim'),[0 0],'color',gray);
title([upper(exp(1)) strrep(exp(2:end),'_',' ') ', high-pass filtered (cutoff=' sprintf('%.2f',cutoff) 'Hz)'],'color','w');
xlabel('Time (TR)');
ylabel(['fMRI signal']);
set(gca,'fontsize',14);
grid on
hold off

for tgi=[1 2];
    
    subplot(3,2,6);
    plot(lags,circularlagcorr(squeeze(nanmean(gdata_lp(1,:,1:9),3))',squeeze(nanmean(gdata_lp(tgi,:,10:18),3))',lags),'linewidth',2);
    hold on
    set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)
    ylim([-0.5 0.5]);
    
    subplot(3,2,5);
    ciplot_claire(squeeze(gdata_lp(tgi,:,:))',keptT,cols(tgi,:),0.3);
    hold on
end
xlim([min(keptT) max(keptT)]);
ylim([-1 1])
gray=[0.7 0.7 0.7];
set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)

line([0 0 ],get(gca,'ylim'),'color',gray);
line(get(gca,'xlim'),[0 0],'color',gray);
title([upper(exp(1)) strrep(exp(2:end),'_',' ') ', low-pass filtered'],'color','w');
xlabel('Time (TR)');
ylabel(['fMRI signal']);
set(gca,'fontsize',14);
grid on
hold off



