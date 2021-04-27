clear all
% close all
set_parameters;

NRFs={'linear', 'log', 'square'};
pauses_sec=[3 1.5 0];
wdN=3000;
scramble=1;
 cols=jet(7);
 
for psi=1%:length(pauses_sec);
    pause_sec=pauses_sec(psi);
    
    for u=3;%2:4
        for vi=0.4;%[0.1:0.3:1.5];
            for speechR=1;%[0.5 1 2]; % 1 is the speech rate of sherlock;
                for af=1%:3;
                    NRF=NRFs{af};
                    load(sprintf('%s/simulation/pause%.1f_m%d_v%.2f_sr%.2f_%s.mat',expdir,pause_sec,u,vi,speechR, NRF),'v_hrf_sec','v','onsetsVectors');
                    
                    [levn,~]=size(v_hrf_sec);
                    
                    figure
                    for lev=1:levn;
                        subplot(2,1,1);
                        plot(v(lev,:),'color',cols(lev+1,:));
                        hold on;
                        xlim([0 size(v,2)])
                        set(gca,'xticklabel',[])
                        ylabel('Neural activity')
                        set(gca,'fontsize',14)
                        subplot(2,1,2);
                        plot(v_hrf_sec(lev,:),'color',cols(lev+1,:));
                        hold on;
                        xlim([0 size(v_hrf_sec,2)])
                        ylabel('BOLD signal')
                        set(gca,'fontsize',14)
                    end
                    xlabel('Time (sec)');
                    legend(cellstr(num2str([1:6]')),'orientation','horizontal')
                    legend boxoff
                    
                    % dow
               
                    % frequency analysis
                    Sxx=[];
                    Syy=[];
                    Sxy=[];
                    for sdi=1:6;
                        for tgi=1:6;
                            [Sxx(sdi,tgi,:) freq] = pwelch(zscore(v_hrf_sec(sdi,:)),50,25,[],1);
                            [Syy(sdi,tgi,:) freq] = pwelch(zscore(v_hrf_sec(tgi,:)),50,25,[],1);
                            [Sxy(sdi,tgi,:) freq] =cpsd(zscore(v_hrf_sec(sdi,:)),zscore(v_hrf_sec(tgi,:)),200,100,[],1);
                            
                        end
                    end
                    coherences=abs(Sxy).^2./(Sxx.*Syy);
                    phases=angle(Sxy);
                    
                    
                    %¡@psd
                    figure;
                    for ni=1:6;
                        plot(log(freq(2:end)),squeeze(Sxx(ni,ni,2:end)),'color',cols(ni+1,:),'linewidth',2);
                        hold on
                    end
                    hold off
                    xlabel('Frequency (Hz)');
                    ylabel('PSD (1/Hz)');
                    xlim(log([freq(2) freq(end)]))
                    set(gca,'xtick',log([0.01 0.04 0.1 0.33]));
                    set(gca,'xticklabels',[0.01 0.04 0.1 0.33]);
                    grid on
                    set(gca,'fontsize',14);
                    title({'Compression ratio:',sprintf('Mean=%d, SD=%.2f',u,sqrt(vi))});
                    legend(cellstr(num2str([1:6]')))
                    legend boxoff
                end
            end
        end
    end
end