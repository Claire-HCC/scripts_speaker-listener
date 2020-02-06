clear all;
close all

set_parameters;
timeUnit='tr' ;
froidir='mor';
lags=-10:-4;
binSize=10;
binStep=1;
rname='HG_L'
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
exp=experiments{2};

fsize=[45 10];
fig1=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
fig2=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);

for ei=1:4;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/LLselfother/regression_LL_binSize' num2str(binSize) '_lag0-0' ],'r2','rnames');
    r2_ll=nanmean(r2,3);
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
    r2_sl=r2;
    
    keptT=[min(find(nansum(r2_sl,1)>0)):max(find(nansum(r2_sl)>0))];

    ri=find(ismember(rnames,rname));
    
    r2_ll_temp=r2_ll(ri,keptT);
    r2_sl_temp=r2_sl(ri,keptT);
    
    set(0, 'currentfigure', fig1);
    subplot(1,4,ei);
    p1=plot([r2_sl_temp'],'r','linewidth',2);
    hold on
    plot([r2_ll_temp'],'b','linewidth',2);
    
    % set(p1.Edge, 'ColorBinding','interpolated', 'ColorData', [uint8(cols(segments_LG,:)*255) uint8(ones(length(keptT),1)*255)].');
    legend({'SL coupling','LL coupling'},'fontsize',13); legend boxoff
    title({exp,rnames{ri}});
    xlim([1 length(keptT)]);
    % ylim([0 1.1*max(max(squeeze(r2_slf(ri,:,:))))]);
    xlabel('Time (TR)'); % ylabel('Coupling');
    set(gca,'fontsize',13);
    
    set(0, 'currentfigure', fig2);
    subplot(1,4,ei);
    scatter(r2_sl_temp',r2_ll_temp',60,'k','filled');
    xlabel('SL coupling'); ylabel('LL coupling');
    [r]=corr(r2_sl_temp',r2_ll_temp')
    title(sprintf('%s, %s, R=%.2f',strrep(rname,'_',' '),exp,r));
    %  if herd_sig_fdr_pos(ri)==1; star='*'; else  star='';end
    %   title(sprintf('Averaged R = %.02f%s',herdm(ri),star));
    set(gca,'fontsize',13);
end



