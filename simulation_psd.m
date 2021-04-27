% the effect of window length
% https://sapienlabs.org/factors-that-impact-power-spectrum-density-estimation/
clear all
% close all
set_parameters;
win=66; %(100/1.5);
crop_start=25;
crop_end=20;
Fs=1/1.5;
simN=30;
tr=1.5;

load(sprintf('%s/simulation/sim%03d/%s',expdir,simN ,'pause0.10_pauseL3.0_m3_v0.50_sr1.00_linear.mat'),'v_hrf_tr');
gdata=v_hrf_tr;
[levn tn]=size(gdata);
  keptT=(crop_start+1):(tn-crop_end);
%% stephens (2013) method
gdata=gdata(:,keptT);

% Stephens: Since the variance is equal across all voxels, all
% spectra have the same integrated area. CHC:sum(Sxx_temp)*(freq(2)-freq(1));  the integrted areas
% would only be exactly the same using periodogram:
gdata=zscore(gdata,0,2);
[Sxx freq] = pwelch(gdata',win,win/2,[],1/tr);
Sxx=Sxx';

save(sprintf('%s/simulation/sim%03d/%s',expdir,simN ,'/psd.mat')  ,'Sxx','keptT','freq','win');

%% plot
fsize=[7 8];
figure('unit','centimeter','position',[0 0 fsize]);
cols=jet(7);
cols=cols([1 3 4 5 6 7],:)*0.95;
fsize=[7 8];
for ni=1:6;
    plot(log(freq(2:end)),Sxx(ni,2:end),'color',cols(ni,:),'linewidth',3);
    %    semilogx(freq(2:end),Sxx(ni,2:end),'color',cols(ni,:),'linewidth',2);
    hold on
end

title('Simulation')
hold off
xlabel('Frequency (Hz)');
ylabel('PSD (1/Hz)');
% xlim(log([min(freq) max(freq)]));
xlim(log([0.01 freq(end)])); % frequency lower than 0.01 (~=1/140) Hz is not meanful after high-pass filtering
ylim([0 35])
%    ylim([0 0.15])
set(gca,'xtick',log([0.01 0.04 0.1 0.33]));
set(gca,'xticklabels',[0.01 0.04 0.1 0.33]);

grid on
set(gca,'fontsize',14)
