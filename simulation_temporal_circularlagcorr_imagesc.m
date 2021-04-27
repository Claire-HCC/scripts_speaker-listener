close all
clear all

% loc='cluster';
set_parameters;
timeUnit='tr' ;
set_parameters;
exp='simulation';
simN=30;

unsig_mask=imread([expdir 'scripts_speaker-listener/unsig_mask.tif']);
dir(sprintf('%s/simulation/sim%03d/*_r.mat',expdir,simN));
u=3;
var=0.5;
speechR=1;
NRF='linear'; % {'linear', 'log','expn', 'triangle','time'};
pauseLen=3;
pauseEffect=0.1;

filelists{1} = cellstr(ls(sprintf('%s/simulation/sim%03d/pause0.10_pauseL3.0_m3_v0.50_sr*_linear_peaks.mat',expdir,simN)));

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
        
        sig=(pfdr_time_peaks<.05 & pfdr_exps_peaks<.05) & (abs(npeaks)<peaks | isnan(npeaks) );
        %  sig=(pfdr_time_peaks<.01 ) & (abs(npeaks)<peaks | isnan(npeaks) );
        peakLags_pfdr=peakLags;
        peakLags_pfdr(sig==0)=NaN;
        disp(filelist{fi})
        peakLags_pfdr
        
        peakLags_bined=discretize(peakLags,[-20 -15 -10 -6 -3 0 1 4  7  11 16 20]);
        % to mask unsignificant peakLags
        im2=imagesc(peakLags_bined,[2 10]);
        colormap jet
        imAlpha=(sig==0)*0.35;
        set(im2,'AlphaData',imAlpha);
        lm=get(gca,'xlim');
        hold on
        h = image([lm(1)-1 lm(2)+1],[lm(1)-1 lm(2)+1],unsig_mask);
        set(h,'AlphaData',max(unsig_mask,[],3)<=85)
        %  uistack(h,'top');
        
        % the significant ones
        im=imagesc(peakLags_bined,[2 10]);
        colormap jet
        imAlpha=ones(size(peakLags));
        imAlpha(sig==0)=0;
        set(im,'AlphaData',imAlpha)
        
        set(gca,'xtick',1:6,'ytick',1:6,'color','k','xlim',lm,'ylim',lm);
        for ni=1:6;
            pos=[ni-0.5 ni-0.5 1 1];
            rectangle('Position',pos ,'edgecolor','w','linewidth',1.5);
        end
        hold off
        
        %         ax=gca;
        %         ax.XAxis.Visible='off';
        %         ax.XAxis.Label.Visible='on';
        %         ax.YAxis.Visible='off';
        %         ax.YAxis.Label.Visible='on';
        %   title(filelist{fi});
        
        yl=ylabel('Seed level');
        %   yl.Position(1)=7.4
        xlabel('Target level');
        set(gca,'fontsize',14)
        
        print(gcf,['temp' num2str(fi) '.tif'],'-dtiff','-r300');
        
    end
end

fs=[];
for fi=1:length(filelist);
    fs=[fs imread(['temp' num2str(fi) '.tif'])];
end
imshow(fs);
imwrite(fs,'temp.tif')