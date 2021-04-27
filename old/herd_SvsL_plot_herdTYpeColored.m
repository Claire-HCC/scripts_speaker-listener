
% clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
lags=-10:-3;
binSize=1;
rname='SOG_L'
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
exp=experiments{2};
% load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_exps_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'b_sl','b_ll','keptT','rnames','lags','herd','herd_p','herd_sig');

for ei=3%1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_leave1out_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
    r2_sl=r2;
    r2_sl_m=nanmean(r2_sl,3);
    
    load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))],'r2','rnames');
    r2_ll=r2;
    r2_ll_m=nanmean(r2,3);
    
    load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_minus1Leave1out_bined/binSize' num2str(binSize) '_lag0-0.mat'],'r2');
    r2_ll0=nanmean(r2,4);
    r2_ll0_m=nanmean(r2_ll0,3);
    
    keptT=[min(find(~isnan(r2_sl(1,:,1)))):max(find(~isnan(r2_sl(1,:,1))))];
    [~,tn,listenerN]=size(r2_sl);
    
    ri=find(ismember(rnames,rname));
    
    r2_ll0_temp=r2_ll0_m(ri,:);
    r2_sl_temp=r2_sl_m(ri,:);
    r2_ll_temp=r2_ll_m(ri,:);
    
    
    herdTypes=nan(tn,1);
    
    for ti=1:length(keptT);
        t=keptT(ti);
        if  r2_sl_temp(t)> r2_ll_temp(t) & r2_ll0_temp(t)>median(r2_ll0_temp(keptT))
            herdTypes(t,1)=1;
        elseif  r2_sl_temp(t)> r2_ll_temp(t) & r2_ll0_temp(t)<=median(r2_ll0_temp(keptT))
            herdTypes(t,1)=2;
        elseif  r2_sl_temp(t)<=r2_ll_temp(t) & r2_ll0_temp(t)<=median(r2_ll0_temp(keptT))
            herdTypes(t,1)=3;
        elseif  r2_sl_temp(t)<= r2_ll_temp(t) & r2_ll0_temp(t)>median(r2_ll0_temp(keptT))
            herdTypes(t,1)=4;
        end
    end
    cols=nan(tn,3);
    col4=[1 0.7 0;0 1 0;0 0 1;0.7 0.1 1];
    cols(keptT,:)=col4(herdTypes(keptT),:);
    
    fsize=[30 9];
    figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
    hold on
    for ti=1:length(keptT);
        
        t=keptT(ti);
        line([t t],[0 1],'color',[cols(t,:) 0.3])
    end
    
    ciplot_claire((squeeze(r2_ll(ri,keptT,:)))',keptT,'r',0.3);
    h = findobj(gca,'Type','line');
    set(h,'linewidth',1.5,'linestyle',':')
    hold on
    ciplot_claire((squeeze(r2_sl(ri,keptT,:)))',keptT,'r',0.3);
    hold on
    ciplot_claire((squeeze(r2_ll0(ri,keptT,:)))',keptT,'b',0.3);
    
    hold off
    title({[upper(exp(1)) exp(2:end) ', ' rnames{ri}]});%,['binSize' num2str(binSize) ', Lags:' num2str(min(lags)) '~' num2str(max(lags)) ]});
    xlim([min(keptT) max(keptT)])
    ylim([0 max(r2_sl_temp(:))+0.05]);
    
    xlabel('Time (TR)');
    ylabel('Neural Coupling');
    set(gca,'fontsize',14);
    h = findobj(gca,'Type','line');
    legend(h,'Listener-Listener','Speaker-Listener','location','northwest');
    legend boxoff
    
    fsize=[9 9];
    figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
    hold on
    for ti=1:length(keptT);
        t=keptT(ti);
        coli=1+ceil(atan2d(r2_sl_temp(t)-median(r2_sl_temp),r2_ll0_temp(t)-median(r2_ll0_temp(keptT))));
        if coli<0; coli=360+coli;end
        scatter(r2_sl_temp(t),r2_ll0_temp(t),20,cols(t,:),'filled');
    end
    xlabel('SL coupling'); ylabel('LL coupling');
    hold off
    title({[upper(exp(1)) exp(2:end) ', ' rnames{ri}]});
    set(gca,'fontsize',14);
    xlabel({'Speaker-Listener coupling','(R-squared)'});
    ylabel({'Listener-Listener coupling','(R-squared)'});
end



