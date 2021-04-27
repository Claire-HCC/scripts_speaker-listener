clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
seeds=networks;
crop_start=25;
crop_end=20;
cutoff=[0.04];
Fs=1/1.5;
filtertypes={'lp','hp'};

for ei=9;%[1 4 11 12 9 10];
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2vox/LL_leave1out/isc'   ],'peakLags','peaks','keptvox');
    maskPeakLag0=zeros(voxn,1);
    thr=sort(peaks(peakLags==0),'descend');
    thr=thr(round(length(keptvox)*0.3));
    maskPeakLag0(keptvox(peaks>thr & peakLags'==0))=1;
    clear peakLags peaks keptvox
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seeds{1} '.mat' ],'gdata','keptvox');
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
    r=nan(length(networks),length(networks),length(lags),length(filtertypes),1);
    
    gdata_lp=[];
    gdata_hp=[];
    for ni=1:6;
        target=networks{ni};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' target '.mat' ],'gdata','keptvox');
        % keep only voxel with peakLag 0
        % disp(sum(ismember(keptvox,find(maskPeakLag0))))
        [~,tn,listenerN]=size(gdata);
        keptT=(crop_start+1):(tn-crop_end);
        
        gdata=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
        gdata(:,:,subjects_excluded{ei})=NaN;
        gdata=nanmean(gdata,1);
        gdata=zscore(gdata(:,keptT,:),0,2);
        gdata_orig(ni,:,:)=gdata;
        for si=1:listenerN;
            gdata_lp(ni,:,si)=filter1('lp',double(gdata(:,:,si)),'fc',cutoff,'fs',Fs);
            gdata_hp(ni,:,si)=filter1('hp',double(gdata(:,:,si)),'fc',cutoff,'fs',Fs);
        end
    end
end



fsize=[30 22];
figure('unit','centimeter','position',[0 0 fsize],'color','k');
nis=[2     4     1     5     3     6];
cols=jet(7);

subplot(3,1,1);
for tgi=[1 6]
    ciplot_claire(squeeze(gdata_orig(nis(tgi),:,:))',keptT,cols(tgi+1,:),0.2);
    
    hold on
end
xlim([min(keptT) max(keptT)]);
ylim([-1.3 1.3])
gray=[0.7 0.7 0.7];
set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)

line([0 0 ],get(gca,'ylim'),'color',gray);
line(get(gca,'xlim'),[0 0],'color',gray);
title([upper(exp(1)) strrep(exp(2:end),'_',' ') ', original'],'color','w');
xlabel('Time (TR)');
ylabel(['fMRI signal']);
set(gca,'fontsize',14);
grid on
hold off

subplot(3,1,2);
for tgi=[1 6]
    ciplot_claire(squeeze(gdata_hp(nis(tgi),:,:))',keptT,cols(tgi+1,:),0.2);
    
    hold on
end
%xlim([50 350]);
ylim([-1 1])
xlim([min(keptT) max(keptT)]);
gray=[0.7 0.7 0.7];
set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)

line([0 0 ],get(gca,'ylim'),'color',gray);
line(get(gca,'xlim'),[0 0],'color',gray);
title([upper(exp(1)) strrep(exp(2:end),'_',' ') ', high-pass filtered'],'color','w');
xlabel('Time (TR)');
ylabel(['fMRI signal']);
set(gca,'fontsize',14);
grid on
hold off

subplot(3,1,3);
for tgi=[1 6];
    ciplot_claire(squeeze(gdata_lp(nis(tgi),:,:))',keptT,cols(tgi+1,:),0.2);
    
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



