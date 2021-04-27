clear all
close all
set_parameters;
timeUnit='tr' ;

froidir='restFc_isc30PercMasked_75Overlap_cluster6';

lags=-15:15;

cols=jet(7)*0.95;
cols=cols([1 3 4 5 6 7],:);
cols(4,:)=[1 0.85 0];
% cols=[ 0.2588    0.5961    0.7098;
%     0.2902    0.1451    0.6667;
%     0.1569    0.2784    0.2039;
%     0.4510    0.3255    0.1137;
%     0.8588    0.0431    0.3569;
%     0.8667    0.5294    0.5529];
gray=[0.9 0.9 0.9];

exps={'sherlock','merlin','pieman_old','black','forgot','ABC','bronx','pieman','pieman_oldWord','pieman_rest','crossExps'};
eis=[1];

for ei=1%:length(eis);
    exp=exps{eis(ei)};
    
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peaks' ],...
        'networks','rzm','rz','rzms','lags','keptT','pfdr_peaks','p_peaks','peakLags','z','peaks','tmid','npeaks');
    data=rzm(:,:,tmid+lags,:);
    
    sig=(pfdr_peaks<(.05) & ~isnan(pfdr_peaks) & ((npeaks)<peaks | isnan(npeaks) ));
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
        
        for ni=1:length(networks);
            temp=squeeze(data(sdi,ni,:));
            peakLag=peakLags(sdi,ni);
            
            if sig(sdi,ni)==1;
                line([peakLag peakLag],[-1 1],'color',cols(ni,:)*0.85)
                % text(peakLag,temp(lags==peakLag),'*','HorizontalAlignment','center','verticalAlignment','middle','fontsize',12,'color',cols(ni,:));
            end
        end
        
        %  set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)
        ylim([-0.1 0.35])
        xlim([min(lags) max(lags)]);
        
        hold off
        
        title([ networks{sdi} ' seed']);%,'color','w')
        set(gca,'xtick',[min(lags) 0 max(lags)],'xticklabel',[min(lags) 0 max(lags)]*exp_parameters.tr(ei));
        if sdi==1;
            xlabel('Lag (sec)');
        end
        set(gca,'fontsize',14)
        grid off
        
        ylabel('R(z)');
        
        print(gcf,['temp' num2str(sdi) '.tif'],'-dtiff','-r1000');
        fs=[fs ;imread(['temp' num2str(sdi) '.tif'])];
        
    end
end
imshow(fs)
imwrite(fs,'temp.tif')