clear all
% close all
set_parameters;
win=66; %(100/1.5);
eis=[1 2 4 11 12 13];
timeUnit='tr';
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
crop_start=25;
crop_end=20;

Fs=1/1.5;
rnames_selected={'vPCUN','HG_L'};
cols=[0 0 1;1 0 0];

eis=[1 2 11 12 ];
fsize=[30 18];
figure('unit','centimeter','position',[0 0 fsize]);
cols=[0 0 1; 1 0 0];
for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    load([expdir '/' exp '/bhv/comprehensionScore.mat'  ],'score');
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/psd.mat'],'Sxx','freq');
    subplot(2,3,eii);
    
    subjs=find(~isnan(squeeze(Sxx(1,1,:))));
    for rii=1:length(rnames_selected);
        
        rname=rnames_selected{rii};
        ri=find(ismember(rnames,rname));
       
        hold on
        scatter(squeeze(sum(Sxx(ri,freq>0.1 & freq<0.2,subjs),2)),score(subjs)',40,cols(rii,:))
        
    end
    
    title([upper(exp(1)) strrep(exp(2:end),'_',' ')],'color','k');
    
    hold off
end

figure;
for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    load([expdir '/' exp '/bhv/comprehensionScore.mat'  ],'score');
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/psd.mat'],'Sxx','freq');
    subplot(2,3,eii);
    
    subjs=find(~isnan(squeeze(Sxx(1,1,:))));
    for rii=1:length(rnames_selected);
        
        rname=rnames_selected{rii};
        ri=find(ismember(rnames,rname));
        for fi=1:length(freq);
            [r(ri,fi) p(ri,fi)]= corr(squeeze(Sxx(ri,fi,subjs)),score(subjs),'type','spearman');
        end
        pfwe(ri,:)=p(ri,:)*length(freq);
        [~,~,pfdr(ri,:)]=fdr(p(ri,:));
        
        plot(freq,r(ri,:),'linewidth',2,'color',cols(rii,:));
        hold on
        scatter(freq(p(ri,:)<0.05),r(ri,p(ri,:)<0.05),40,cols(rii,:));
    end
    
end


