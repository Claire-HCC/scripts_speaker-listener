close all
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='restFc_isc30PercMasked_75Overlap_cluster6_audPausesResid';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};

lags=-15:15;
exps={'sherlock','merlin','pieman_old','black','forgot','ABC','bronx','pieman','pieman_oldWord','pieman_rest','crossExps'};
eis=[1:11];
fs=[];

for ei=1%:length(eis);
    exp=exps{eis(ei)};
    
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peaks' ],...
        'keptT','networks','tmid','rzm','rzm_timeReversed','peakLags','pfdr_peaks','peaks','npeaks');
    sig=(pfdr_peaks<(.05) & ~isnan(pfdr_peaks) & (abs(npeaks)<peaks | isnan(npeaks) ));
    
    
    for sdi=1:length(networks);
        
        fsize=[15, 4.5];
        h=figure('unit','centimeter','position',[0 0 fsize]);
        
        %   subplot(length(networks),1,sdi);
        [~,~,tn]=size(rzm);
        
        cols=gray(100);
        cols=cols(discretize(log2([1:0.01:3.5]+1),100),:);
        colormap(cols)
        imagesc(zscore(squeeze(rzm(sdi,:,:)),0,2),[1 3.5]);
        
        xlim([tmid+[min(lags)-0.5 max(lags)+0.5]]);
        set(gca,'xtick',tmid+[min(lags) 0 max(lags)],'xticklabels',strrep(cellstr(num2str(1.5*[min(lags) 0 max(lags)]')),' ',''));
        set(gca,'ytick',1:length(networks),'yticklabels',networks);
        
        cols=jet(64)*0.95;
        
        for tgi=1:length(networks);
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
        ylabel([networks{sdi} ' seed'],'FontWeight','bold','fontsize',16)
        xlabel('Lag (sec)','FontWeight','bold','fontsize',16,'fontweight','bold');
        
        print(gcf,['temp' num2str(sdi) '.png'],'-dpng','-r1000');
        fs=[fs ;imread(['temp' num2str(sdi) '.png'])];
    end
    
    % c = colorbar;
    % c.Label.String = {'R (z)','normalized across lags'};
    %  title({strrep([upper(exp(1)) exp(2:end)],'_',' ')});
end
imshow(fs);
imwrite(fs,'temp.tif')

