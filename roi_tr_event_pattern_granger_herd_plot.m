% bined r2 vs std between subjects
clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
rname='vPCUN';
% -4 for merlin, -3 for sherlock
lags=-60:-4;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

for ei=3%:4;%1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_event_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2_l','r2_ls','r2_s','eventLabels_LG_kept','rnames');
    r2_s2l=r2_ls-r2_l;
    r2_l2l=r2_ls-r2_s;
    eventLabels_LG_kept_sl=eventLabels_LG_kept;
    
  load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_event_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm'],'p_s2l');

    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_event_lag0-0' ],'r2','eventLabels_LG_kept');
    r2_ll=nanmean(r2,4);
    eventLabels_LG_kept_ll=eventLabels_LG_kept;
    
    % try dividing the story into early and late parts
    eventLabels_LG_kept=intersect(eventLabels_LG_kept_sl,eventLabels_LG_kept_ll);
    
    evs=(unique(eventLabels_LG_kept));
    % cols=distinguishable_colors(length(evs));
    cols=hot(length(evs)+10);

    % ris=find(ismember(rnames,{'HG_L','vPCUN','aCUN'}));
    ris=find(ismember(rnames,rname));
    % [~,ri]=max(herdm);
    
    for rii=1:length(ris);
        ri=ris(rii);
        
        fsize=[20 9];
        figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
        
        r2_ll_g=nanmean(r2_ll,3);
        
        subplot(1,2,1);
        hold on
        p1=plot([r2_s2l(ri,evs)'],'r','linewidth',2);
        plot([r2_ll_g(ri,evs)'],'k','linewidth',2);
       % area(p_s2l(ri,evs)==0)
        hold off
        legend({'SL coupling','LL coupling'},'fontsize',13); legend boxoff
        title({exp,rnames{ri}});
        xlim([1 length(evs)]);
        % ylim([0 1.1*max(max(squeeze(r2_slf(ri,:,:))))]);
        xlabel('event'); % ylabel('Coupling');
        set(gca,'fontsize',13);
        
        subplot(1,2,2);
      %  scatter(r2_s2l(ri,evs)',r2_ll_g(ri,evs)',40,cols(evs,:),'filled');
       scatter(r2_s2l(ri,evs)',r2_ll_g(ri,evs)',40,'k','filled');
        hold on
        evs_sig=intersect(evs,find(p_s2l(ri,:)==0));
        scatter(r2_s2l(ri,evs_sig)',r2_ll_g(ri,evs_sig)',40,'r','filled');
        
        [r p]=corr(r2_s2l(ri,evs)',r2_ll_g(ri,evs)');
        xlabel('SL coupling'); ylabel('LL coupling');
        
        % if herd_sig_fdr_pos(ri)==1; star='*'; else  star='';end
        % title(sprintf('Averaged R = %.02f%s',herdm(ri),star));
        
        title(sprintf('Averaged R = %.2f, p=%.3f',r,p));
        set(gca,'fontsize',13);
    end
    
    clear herd herd_null herd_p herd_sig_fdr_pos herd_sig_fdr_neg
end
