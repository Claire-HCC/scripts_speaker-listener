close all
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;

simN=30;
% 
cols=jet(7)*0.95;
cols=cols([1 3 4 5 6 7],:);
cols(4,:)=[1 0.85 0];
% cols=[ 0.2588    0.5961    0.7098;
%     0.2902    0.1451    0.6667;
%     0.1569    0.2784    0.2039;
%     0.4510    0.3255    0.1137;
%     0.8588    0.0431    0.3569;
%     0.8667    0.5294    0.5529];

NRFs={'linear', 'log', 'triangle','time'};
networks=cellstr(num2str([1:6]'));
tr=1.5;
gray=[0.7 0.7 0.7];
lags=-15:15;

filelists{1} = cellstr(ls(sprintf('%s/simulation/sim%03d/pause0.10_pauseL3.0_m3_v0.50_sr1.00_linear_peaks.mat',expdir,simN)));

% filelists{1} = cellstr(ls(sprintf('%s/simulation/sim%03d/pause*_pauseL3.0_m3_v0.50_sr1.00_linear_peaks.mat',expdir,simN)));
% filelists{2} = cellstr(ls(sprintf('%s/simulation/sim%03d/pause0.10_pauseL*_m3_v0.50_sr1.00_linear_peaks.mat',expdir,simN)));
% filelists{3} = cellstr(ls(sprintf('%s/simulation/sim%03d/pause0.10_pauseL3.0_m*_v0.50_sr1.00_linear_peaks.mat',expdir,simN)));
% filelists{4} = cellstr(ls(sprintf('%s/simulation/sim%03d/pause0.10_pauseL3.0_m3_v*_sr1.00_linear_peaks.mat',expdir,simN)));
% filelists{5} = cellstr(ls(sprintf('%s/simulation/sim%03d/pause0.10_pauseL3.0_m3_v0.50_sr*_linear_peaks.mat',expdir,simN)));
% filelists{6} = cellstr(ls(sprintf('%s/simulation/sim%03d/pause0.10_pauseL3.0_m3_v0.50_sr1.00_*_peaks.mat',expdir,simN)));

for parami=1:length(filelists);
    filelist=filelists{parami};
    fsize=[ 7 7.3];
    figure('unit','centimeter','position',[0 0 fsize],'InvertHardcopy', 'off','color','w');
    
    fs=[];
    for fi=1:length(filelist);
        load(sprintf('%s/simulation/sim%03d/%s',expdir,simN ,filelist{fi}));
        data=rzm(:,:,tmid+lags,:);
        
        sig=(pfdr_time_peaks<.05 & pfdr_exps_peaks<.05) & (abs(npeaks)<peaks | isnan(npeaks) );
        fs=[];
        for sdi=1:length(networks);
            seed=networks{sdi};
            fsize=[7 7];
            figure('unit','centimeter','position',[0 0 fsize]);%,'color','k','InvertHardcopy', 'off');%
            
            line([min(lags) max(lags)],[0 0],'color',gray)
            hold on;
            for ni=1:length(networks);
                temp=squeeze(data(sdi,ni,:,:));
                %     ciplot_claire(temp',lags,cols(ni,:),0.3);
                plot(lags,tanh(nanmean(temp,2)),'color',cols(ni,:),'linewidth',2);
                hold on
            end
            
               ylim([-0.4 1.2])
            xlim([min(lags) max(lags)]);
            
            for ni=1:length(networks);
                temp=squeeze(data(sdi,ni,:));
                peakLag=peakLags(sdi,ni);
                
                if sig(sdi,ni)==1;
                    line([peakLag peakLag],get(gca,'ylim'),'color',cols(ni,:)*0.85)
                    % text(peakLag,temp(lags==peakLag),'*','HorizontalAlignment','center','verticalAlignment','middle','fontsize',12,'color',cols(ni,:));
                end
            end
            
            %  set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)
         
            
            hold off
            
            title([ 'Level ' num2str(sdi) ' seed']);%,'color','w')
            set(gca,'xtick',[min(lags) 0 max(lags)],'xticklabel',[min(lags) 0 max(lags)]*tr);
            if sdi==1;
                xlabel('Lag (sec)');
            end
            set(gca,'fontsize',14)
            grid off
            
            ylabel('R');
            
            print(gcf,['temp' num2str(sdi) '.tif'],'-dtiff','-r100');
            fs=[fs ;imread(['temp' num2str(sdi) '.tif'])];
            
        end
    end
    imshow(fs)
    imwrite(fs,'temp.tif')
end