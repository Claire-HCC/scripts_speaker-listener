% the effect of window length
% https://sapienlabs.org/factors-that-impact-power-spectrum-density-estimation/
clear all
close all
set_parameters;
win=66; %(100/1.5);
eis=[1 2 4 11 12];
for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    subplot(2,3,eii);
    load([expdir exp '\fmri\timeseries\tr\roi\mor\listenerAll_vPCUN.mat'])
    g=nanmean(zscore(gdata,0,2),3);
    [pxx, f] = pwelch(g',win,win/2,[],1/tr(ei));
    plot(f,nanmean(pxx,2),'linewidth',2);
    hold on
    load([expdir exp '\fmri\timeseries\tr\roi\mor\listenerAll_HG_L.mat'])
    g=nanmean(zscore(gdata,0,2),3);
    [pxx, f] = pwelch(g',win,win/2,[],1/tr(ei));
    plot(f,nanmean(pxx,2),'linewidth',2);
    grid on;
    title(exp)
    hold off
end
legend('vPCUN','HG L')