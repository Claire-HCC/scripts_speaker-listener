clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;

roi='isc_peak_precuneus';
crop_start=25;
crop_end=20;
perc=0.15;
lags=-20:20;
cols=jet(length(lags));

for ei=[5  ];%[1 4 11 12 9 10];
    exp=exp_parameters.experiments{ei};
    
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/LL_leave1out/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peaks' ],...
        'keptvox','peakLags','pfdr_peaks');
    keptvox_mask=keptvox(pfdr_peaks<.05);
    peakLags=peakLags(pfdr_peaks<.05);
  %  peakLags_uni=unique(peakLags);
    
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll_isc30PercMasked.mat',expdir,exp,timeUnit);
    load(f,'gdata','keptvox');
    gdata(:,:,exp_parameters.subjects_excluded{ei})=NaN;
    
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    
     peakLags_uni=unique(peakLags);
    peakLags_uni=[-13 0];
    figure;
    for pi=1:length(peakLags_uni);
        peakLag=peakLags_uni(pi);
        
        vis=ismember(keptvox,keptvox_mask(peakLags==peakLag));
        x=nanmean(gdata(vis,:,:),1);
        x=zscore(x(:,keptT,:),0,2);
        
        %     ciplot_claire(squeeze(gdata_orig(ni,:,:))',keptT,cols(find(ismember(network_newOrder,network))+1,:),0.2);
        ciplot_claire(squeeze(x)',keptT,cols(ismember(lags,peakLag),:),0);
        hold on;
    end
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
