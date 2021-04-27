% the effect of window length
% https://sapienlabs.org/factors-that-impact-power-spectrum-density-estimation/
clear all
close all
set_parameters;


timeUnit='tr';
crop_start=25;
crop_end=20;
froidir='restFc_isc30PercMasked_75Overlap_cluster6';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};
Fs=1/1.5;

%% plot
cols=jet(7);
cols=cols([1 3 4 5 6 7],:)*0.95;
fsize=[30 27];
figure('unit','centimeter','position',[0 0 fsize]);
Sxxs=[];
for ei=1:10;
    exp=exp_parameters.experiments{ei};
    
    load([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/network/' froidir '/L_g/networks_psd'  ],'Sxx','networks','keptT','freq');
    nansum(Sxx,2)'
    
    subplot(3,4,ei);
    if ei==1;
        Sxxs=Sxx;
    else
        Sxxs(:,:,ei)=Sxx;
    end
    
    for ni=1:6;
        plot(log(freq(2:end)),Sxx(ni,2:end),'color',cols(ni,:),'linewidth',2);
        %    semilogx(freq(2:end),Sxx(ni,2:end),'color',cols(ni,:),'linewidth',2);
        hold on
    end
    
    title([upper(exp(1)) strrep(exp(2:end),'_',' ')])
    hold off
    xlabel('Frequency (Hz)');
    ylabel('PSD (1/Hz)');
    xlim(log([min(freq) max(freq)]));
    %  xlim(log([0.01 freq(end)])); % frequency lower than 0.01 (~=1/140) Hz is not meanful after high-pass filtering
    ylim([0 35])
    %    ylim([0 0.15])
    set(gca,'xtick',log([0.01 0.04 0.1 0.33]));
    set(gca,'xticklabels',[0.01 0.04 0.1 0.33]);
    
    grid on
    set(gca,'fontsize',14)
end


% subplot(3,4,ei+1);
fsize=[7 8];
figure('unit','centimeter','position',[0 0 fsize]);
Sxx=mean(Sxxs,3);
for ni=1:6;
    plot(log(freq(2:end)),Sxx(ni,2:end),'color',cols(ni,:),'linewidth',3);
    %    semilogx(freq(2:end),Sxx(ni,2:end),'color',cols(ni,:),'linewidth',2);
    hold on
end

title('Story (N=8)')
hold off
xlabel('Frequency (Hz)');
ylabel('PSD (1/Hz)');
% xlim(log([min(freq) max(freq)]));
xlim(log([0.01 freq(end)])); % frequency lower than 0.01 (~=1/140) Hz is not meanful after high-pass filtering
ylim([0 25])
%    ylim([0 0.15])
set(gca,'xtick',log([0.01 0.04 0.1 0.33]));
set(gca,'xticklabels',[0.01 0.04 0.1 0.33]);

grid on
set(gca,'fontsize',14)
