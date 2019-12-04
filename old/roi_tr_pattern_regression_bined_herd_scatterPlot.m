
clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
lags=-4:4;

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
exp=experiments{2};
% load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_exps_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'b_sl','b_ll','keptT','rnames','lags','herd','herd_p','herd_sig');

ris=[3 54];
% ris=find(sum(herd_sig==0,2)==0);
for i=1:length(ris);
    ri=ris(i);
    
    fsize=[35 20];
    fig1=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
    fig2=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
    for ei=1:4;
        exp=experiments{ei};
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_bined'],'couplingz','keptT','rnames');
        cp_ll=nanmean(couplingz,3);
        
        
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_bined'],'couplingz','keptT');
        cp_sl=couplingz;
        
        
        rname=rnames{ri};
        
        figure(fig1);
        subplot(2,2,ei);
        [r p ]=corr(cp_sl(ri,:)',cp_ll(ri,:)');
        
        p = polyfit(cp_sl(ri,:)',cp_ll(ri,:)', 1);
        Y = polyval(p, (cp_sl(ri,:)'));
        
        scatter(cp_sl(ri,:),cp_ll(ri,:),20,'k','filled');
        hold on
        plot(cp_sl(ri,:),Y,'color',[0.5 0.5 0.5],'linewidth',1.5);
        title([rname ', ' exp ', ' sprintf('R=%.02f',r)]);
        hold off
        xlabel('Speaker-Listerner Coupling');
        ylabel('Listener-Listerner Coupling');
        
        figure(fig2);
        subplot(2,2,ei);
        plot(zscore([cp_sl(ri,:)' cp_ll(ri,:)'],0,1));
        legend({'Speaker-Listerner Coupling','Listener-Listerner Coupling'});
        legend boxoff
        xlim([0 length(cp_sl(ri,:))]);
        title([rname ', ' exp ', ' sprintf('R=%.02f',r)]);
        xlabel('Time (TR)');
        ylabel('zscored coupling')
        
    end
    
end


