%% cross stories LL bined consistency, networks2networks
clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='restFc_isc15PercMasked_50Overlap_cluster5';

lags=-20:20;

eis=1:9;
for i=1:length(eis);
    ei=eis(i);
    exp=exp_parameters.experiments{ei};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peaks' ],...
        'networks','rz','rzm','lags','keptT','pfdr_peaks','p_peaks','peakLags','peaks','pfdr_npeaks');
    
    peakLags_ll(:,:,ei)=peakLags;
    peaks_ll(:,:,ei)=peaks;
    
    sigs(:,:,ei)=(pfdr_peaks<(.05) & ~isnan(p_peaks) & (pfdr_npeaks>.05 | isnan(pfdr_npeaks) ));
    
end

fsize=[40 8 ];
figure('unit','centimeter','position',[0 0 fsize]);
subploti=reshape(1:(length(eis)*(length(eis)-1)),length(eis),length(eis)-1)';


for sdi=1:length(networks);
    
    subplot(1,length(networks),sdi);
    
    x=repmat([1:(length(networks))]',1,size(peakLags_ll,3));
    y=squeeze(peakLags_ll(sdi,:,:)).*1.5;
    sig=squeeze(sigs(sdi,:,:));
    
    
    t = 2*pi*rand(size(x));
    r = 0.3*sqrt(rand(size(x)));
    x_jittered = x + r.*cos(t);
    y_jittered = y + r.*sin(t);
    
    y_log=sign(y).*log(abs(y)+1);
    
    [r1 p1]=corr(x(sig),y(sig ),'type','spearman','tail','right');
    scatter(x_jittered(sig(:,1:8)),y_log(sig(:,1:8)),15,[0.8 0.1 0.1],'filled','MarkerFaceAlpha',1,'MarkerEdgeColor',[1 0 0]);
    
    hold on
    % [r2 p2]=corr(x(shared & triui),y(shared & triui),'type','spearman','tail','right');
    
    scatter(x_jittered(sig(:,1:8)==0),y_log(sig(:,1:8)==0),10,'k','filled','MarkerFaceAlpha',0.2);
    lags_log=sign(lags*1.5).*log(abs(lags*1.5)+1);
   
    set(gca,'xtick',[1:(length(networks))],'xticklabel',networks);
    xtickangle(45)
   
    set(gca,'ytick',lags_log(1:5:end),'yticklabel',lags(1:5:end)*1.5);
   
    xlim([0 length(networks)]);
    ylim([min(lags_log) max(lags_log)]);
    set(gca,'fontsize',12);
    title({[networks{sdi} ' seed'],sprintf('R=%.2f (N=%d)',r1,sum(sig(:)))})
      
 
    xlabel('Target networks');
    ylabel('Peak lag (sec)')
    % xlabel({[upper(exp1(1)) strrep(exp1(2:end),'_',' ') ]},'fontsize',12,'fontweight','bold');
    % ylabel({[upper(exp2(1)) strrep(exp2(2:end),'_',' ') ]},'fontsize',12,'fontweight','bold');
end

