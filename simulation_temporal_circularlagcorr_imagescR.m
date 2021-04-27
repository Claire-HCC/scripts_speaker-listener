close all
clear all

% loc='cluster';
set_parameters;
timeUnit='tr' ;
set_parameters;
exp='simulation';
simN=30;
lags=-15:15;
% u=3;
% var=0.5;
% speechR=0.5;
% NRF='linear'; % {'linear', 'log','expn', 'triangle','time'};
% pauseLen=3;
% pauseEffect=0.1;

filelists{1} = cellstr(ls(sprintf('%s/simulation/sim%03d/pause0.10_pauseL3.0_m3_v0.50_sr1.00_linear_peaks.mat',expdir,simN)));

for parami=1:length(filelists);
    filelist=filelists{parami};
    
    for fi=1:length(filelist);
        load(sprintf('%s/simulation/sim%03d/%s',expdir,simN ,filelist{fi}),'rzm','tmid','peakLags','npeaks','peaks','pfdr_exps_peaks','pfdr_time_peaks');
        sig=(pfdr_time_peaks<.05 & pfdr_exps_peaks<.05) & (abs(npeaks)<peaks | isnan(npeaks) );
        [levn,~,tn,~]=size(rzm);
        
        fs=[];
        for sdi=1 :levn;
            
            fsize=[15, 4.5];
            h=figure('unit','centimeter','position',[0 0 fsize]);
            
            %   subplot(length(networks),1,sdi);
            [~,~,tn]=size(rzm);
            cols=gray(100);
            
            cols=cols(discretize(log2([1:0.01:5]+1),100),:);
            colormap(cols)
            imagesc(zscore(squeeze(rzm(sdi,:,:)),0,2),[1 5]);
            
            xlim([tmid+[min(lags)-0.5 max(lags)+0.5]]);
            
            set(gca,'xtick',tmid+[-15 0 15],'xticklabels',{'-22.5','0','22.5'});
            set(gca,'ytick',1:levn,'yticklabels',cellfun(@(x) ['Level ' x],cellstr(num2str([1:6]')),'Uniformoutput',0));
            
            cols=jet(64)*0.95;
            
            for tgi=1:levn;
                pos=[tmid+peakLags(sdi,tgi) tgi-0.5 0 1];
                
                if sig(sdi,tgi)==1;
                    peakLag=peakLags(sdi,tgi);
                    coli=round(discretize(peakLag,[ -15 -10 -6 -3 0 1 4  7  11 16])/9*64);
                    line([peakLag peakLag]+tmid,[tgi-0.5 tgi+0.5],'color',cols(coli,:),'linewidth',8)
                    
                    
                elseif sig(sdi,tgi)==0 & ~isnan(peakLags(sdi,tgi));
                    peakLag=peakLags(sdi,tgi);
                    coli=round(discretize(peakLag,[ -15 -10 -6 -3 0 1 4  7  11 16])/9*64);
                    line([peakLag peakLag]+tmid,[tgi-0.5 tgi+0.5],'color',cols(coli,:),'linewidth',3)
                    
                end
                
            end
            
            hold off
            set(gca,'fontsize',12)
            ylabel([ 'Level ' num2str(sdi) ' seed'],'FontWeight','bold','fontsize',16)
            xlabel('Lag (sec)','FontWeight','bold','fontsize',16,'fontweight','bold');
            
            print(gcf,['temp' num2str(sdi) '.png'],'-dpng','-r1000');
            fs=[fs ;imread(['temp' num2str(sdi) '.png'])];
        end
        
    end
end

imshow(fs)
imwrite(fs,'temp.tif')
