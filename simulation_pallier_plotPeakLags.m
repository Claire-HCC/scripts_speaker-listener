clear all
close all
set_parameters;

NRFs={'linear', 'log', 'square'};
pauses_sec=[3 1.5 0];
wdN=3000;
scramble=1;
simN=1;
cols=jet(7);


for psi=1;%1:length(pauses_sec);
    pause_sec=pauses_sec(psi);
    
    for speechR=1;%[0.5 1 2]; % 1 is the speech rate of sherlock;
        
        for af=1%:3;
            NRF=NRFs{af};
            fsize=[27 20];
            figure('unit','centimeter','position',[0 0 fsize]);
            spi=1;
            
            for u=2:4
                for vi=[0.1:0.3:1.5];%[0.1:0.3:1.5];
                    
                    
                    load(sprintf('%s/simulation/sim%03d/pause%.1f_m%d_v%.2f_sr%.2f_%s_stats.mat',expdir,simN,pause_sec,u,vi,speechR, NRF),'peaks','peakLags','peaks_pfwe','peakLags_pfwe','rzm','');
                    
                    levn=size(peaks,1);
                    % peak lag imagesc
                    
                    subplot(3,5,spi);
                    
                    %            peakLags_bined=discretize(peakLags,[-20 -15 -10 -6 -3 0 1 4  7  11 16 20]*1.5);
                    %    im=imagesc(peakLags_bined,[2 10]);
                    %   peakLags_pfwe
                    %  im=imagesc(peakLags,[-5 5]);
                    %  set(gca,'color','k')
                    % set(im,'AlphaData',~isnan(peakLags_pfwe));
                    
                    im=imagesc(squeeze(rzm(:,:,lags==0)),[-1 1]);
                    colormap jet
                    
                    title({'Compression ratio:',sprintf('Mean=%d, SD=%.2f',u,sqrt(vi))});
                    
                    spi=spi+1;
                end
            end
        end
    end
end