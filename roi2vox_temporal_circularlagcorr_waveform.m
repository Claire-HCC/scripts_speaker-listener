clear all
% close all

set_parameters
timeUnit='tr' ;
froidir='isc_peak';
lags=-20:20;
rois={'precuneus','Angular_L'};
perc=0.3;
cols=jet(length(lags));
gray=[0.9 0.9 0.9];

for ei=5;%[1:10];
    exp=exp_parameters.experiments{ei};
    
    for ri=1;%1:length(rois);
        roi=rois{ri};
        load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_leave1out/' roi '2voxIsc' num2str(perc*100) 'PercMasked_audResid_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peaks' ],...
            'keptvox','peakLags','pfdr_peaks','rzm','pfdr','pfdr_npeaks');
    %    sig=(pfdr_peaks<.05 & (pfdr_npeaks>.05 | isnan(pfdr_npeaks) ));
     sig=(pfdr_npeaks<.05 );
        
        [~,tn]=size(rzm);
        tmid=(tn-1)/2+1;
        
        figure;
        [~,vis]=sort(peakLags);
        vis=vis(sig(vis)==1);
        imagesc(rzm(vis,:));
        
        peakLags_uni=unique(peakLags(sig));
        %  peakLags_uni=[-13 0];
        % figure('unit','centimeter','position',[0 0 fsize],'color','k');
        figure;
        for pi=1:length(peakLags_uni);
            peakLag=peakLags_uni(pi);
            vis=(peakLags==peakLag & sig);
            
            % subplot(5,5,find(lags==peakLag));
            %  imagesc(rzm(vis,:),[-0.2 0.2]);
            
            line([min(lags) max(lags)],[0 0],'color',gray)
            hold on;
            
            temp=squeeze(nanmean(nanmean(rzm(vis,tmid+lags,:),1),3));
            plot(lags,temp,'color',[cols(lags==peakLag,:) 0.8],'linewidth',1.5);
            hold on
        end
        %   set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)
        ylim([-0.2 0.3])
        xlim([min(lags) max(lags)]);
        
        hold off
        grid on
        set(gca,'xtick',[-20:10:20],'xticklabel',[-20:10:20]*exp_parameters.tr(ei));
        xlabel('Lag (sec)');
        set(gca,'fontsize',10)
        
        if peakLag==min(lags);
            ylabel({[upper(exp(1)) strrep(exp(2:end),'_',' ')],'R(z)'})
        else
            ylabel('R(z)');
        end
        title([upper(exp(1)) exp(2:end)]);
    end
end


