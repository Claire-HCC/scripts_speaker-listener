close all
clear all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
lags=-15:15;
type='env';
role='listener';
crop_start=0;
crop_end=0;
eis=[1 2 4 11 12];% 1:2;%1:4;

fsize=[40 7];
figure('unit','centimeter','position',[0 0 fsize]);
for i=1;%1:length(eis);
    ei=eis(i);
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2vox/LL_leave1out/isc'   ],'peakLags','peaks','keptvox');
    maskPeakLag0=zeros(voxn,1);
    thr=sort(peaks(peakLags==0),'descend');
    thr=thr(round(length(keptvox)*0.3));
    maskPeakLag0(keptvox(peaks>thr & peakLags'==0))=1;
    clear peakLags peaks keptvox
    
    subplot(1,5,i)
    f=ls([expdir exp '/sound/*_listener_audenv.mat' ]);
    load([expdir experiments{ei} '/sound/' f ]);
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_HG_L.mat' ],'gdata','keptvox');
   %  load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_Auditory_Language.mat' ],'gdata','keptvox');
    gdata=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
    gdata(:,:,subjects_excluded{ei})=NaN;
    gdata=nanmean(gdata,1);
    
    [~,tn,listenerN]=size(gdata);
    tn=min(min(tn,length(aud)));
    
    r1=(lagcorr_claire(aud(1:(tn-crop_start)),nanmean(gdata((crop_start+1):tn),3)',lags))';
    
    plot(lags*tr(ei),r1,'linewidth',2,'color','k');
    xlabel('Lag (sec)');
    ylabel('R');
    title([ upper(exp(1)) strrep(exp(2:end),'_',' ')],'fontsize',14);
    set(gca,'fontsize',14);
    
    xlim([min(lags) max(lags)])
 %   ylim([-0.2 0.2])
 %   set(gca,'xticklabels',[-20:5:20]*tr(ei),'xtick',-20:4:20);
    grid on
end