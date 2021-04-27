% the effect of window length
% https://sapienlabs.org/factors-that-impact-power-spectrum-density-estimation/
clear all
close all
set_parameters;
win=100; %(100/1.5);
eis=[1 2 4 11 12 9 10];
froidir='mor';
crop_start=25;
crop_end=20;


for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    
    load([ expdir exp '\sound\' exp '_listener_audhrf.mat'],'aud');
    
    load([expdir exp '\fmri\timeseries\tr\network\'  froidir '\listenerAll_DMN1.mat']);
    
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    aud=zscore(aud(keptT));
    [Sxx, freq] = pwelch(aud,win,win/2,[],1/tr(ei));
    
    save([expdir '/' exp '/sound/audhrf_psd'  ],'Sxx','keptT','freq');
    clear Sxx
end


eis=[1 2 4 11 12  ];
fsize=[30 18];
figure('unit','centimeter','position',[0 0 fsize]);
for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    load([expdir '/' exp '/sound/audhrf_psd'  ],'Sxx','networks','keptT','freq');
    
    subplot(2,3,eii)
    
    plot(freq,Sxx,'linewidth',2,'color','k')
    hold on
    
    grid on;
    xlim([0 0.1]);
    title([upper(exp(1)) strrep(exp(2:end),'_',' ')])
    hold off
    
    clear pxx
end