% the effect of window length
% https://sapienlabs.org/factors-that-impact-power-spectrum-density-estimation/
clear all
% close all
% loc='cluster';
set_parameters;
win=66; %(100/1.5);
timeUnit='tr';
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
crop_start=25;
crop_end=20;
Fs=1/1.5;
network_newOrder={
    'Auditory_Language',...
    'DMN2',...
    'Attention',...
    'Executive',...
    'DMN1',...
    'Visual'...
    };
cols=jet(7);

eis=[1 2 4 11 12 13 9 10];
cutoff=0.04;
fsize=[30 27];
figure('unit','centimeter','position',[0 0 fsize]);

for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/roi/' froidir '/LL_leave1out/rois2rois_cpsd'  ],'rnames','keptT','freq','win','Sxx','Syy','Sxy');
    
    ris=[];
    for ni=1:length(network_newOrder);
        ris=[ris ; find(ismember(rnames,roi_table.region(ismember(roi_table.network,network_newOrder{ni}))))];
    end
    
    subplot(3,3,eii);
    coherences=abs(Sxy).^2./(Sxx.*Syy);
    
    for ri=1:length(rnames);
        if ismember(ri,[1 11]);
            rname=rnames{ri};
            network=roi_table.network(ri);
            ni=find(ismember(network_newOrder,network));
            c1=squeeze(nanmean(coherences(ri,ri,freq>0 & freq<=cutoff,:)));
            c2=squeeze(nanmean(coherences(ri,ri,freq>cutoff,:)));
            
            c1=c1(~isnan(c1));
            c2=c2(~isnan(c2));
           
           errorbar(nanmean([c1 c2]),ci([c1 c2],0.95)/2,'linewidth',2,'color',cols(ni+1,:));
            hold on
        end
    end
    hold off
    xlim([0.5 2.5])
    
    set(gca,'xtick',1:2,'xticklabels',{'< 0.04 Hz','> 0.04 Hz'});
    
    ylabel({'Between-subject coherence'})
    xlabel('Frequency')
    
    set(gca,'fontsize',14)
    title([upper(exp(1)) strrep(exp(2:end),'_',' ')],'color','k');
    grid on
    
end
