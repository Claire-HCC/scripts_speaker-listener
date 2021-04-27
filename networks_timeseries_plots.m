clear all
close all

set_parameters;
timeUnit='tr' ;

froidir='restFc_isc30PercMasked_75Overlap_cluster6';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};
% cols=jet(7);
% cols=cols([1 3 4 5 6 7],:);
% cols(4,:)=[1 0.85 0];


crop_start=25;
crop_end=20;
seed='AUD';
for ei=[1 ];%[1 4 11 12 9 10];
    exp=exp_parameters.experiments{ei};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed ],'gdata','keptvox');
    
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    
    gdata(:,:,exp_parameters.subjects_excluded{ei})=NaN;
    gdata=nanmean(gdata,1);
    gdata_seed=zscore(gdata(:,keptT,:),0,2);
    
    figure('unit','centimeter','position',[0 0 6 30]);
    
    spi=0;
    for ni=1:6;%length(networks);
        spi=spi+1;
        subplot(6,1,spi);
        %  network=networks{ni};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' networks{ni} ],'gdata','keptvox');
        
        [~,tn,listenerN]=size(gdata);
        keptT=(crop_start+1):(tn-crop_end);
        
        gdata(:,:,exp_parameters.subjects_excluded{ei})=NaN;
        gdata=nanmean(gdata,1);
        gdata=zscore(gdata(:,keptT,:),0,2);
        plot(keptT*exp_parameters.tr(ei),nanmean(squeeze(gdata_seed(:,:,1:ceil(listenerN/2)))'),'color','k','linewidth',4);
        hold on
        plot(keptT*exp_parameters.tr(ei),nanmean(squeeze(gdata(:,:,(ceil(listenerN/2)+1):end))'),'color',[0 0 0 0.5],'linewidth',4);
        hold on;
        
        xlabel('Time (sec)');
        ylabel({'fMRI signal'});
        set(gca,'fontsize',14);
        hold off
    %    set(gca,'xtick',[],'ytick',[])
       xlim([400 430]);
        legend({[seed ' seed'],[networks{ni} ' target']},'fontsize',12,'location','northwest');
        legend boxoff
        ylim([-0.4 1.1]);
        %  xlim([500 520]);
    end
end

