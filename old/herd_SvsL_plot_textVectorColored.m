
% clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
lags=-10:-3;
binSize=1;
rname='HG_R'
kws={'na','arousal'};

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
exp=experiments{2};
% load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_exps_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'b_sl','b_ll','keptT','rnames','lags','herd','herd_p','herd_sig');
text_trVectorShift=1;
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
    
%     herdTypes=nan(tn,1);
%     
%     for ti=1:length(keptT);
%         t=keptT(ti);
%         if  r2_sl_temp(t)> r2_ll_temp(t) & r2_ll0_temp(t)>median(r2_ll0_temp(keptT))
%             herdTypes(t,1)=1;
%         elseif  r2_sl_temp(t)> r2_ll_temp(t) & r2_ll0_temp(t)<=median(r2_ll0_temp(keptT))
%             herdTypes(t,1)=2;
%         elseif  r2_sl_temp(t)<=r2_ll_temp(t) & r2_ll0_temp(t)<=median(r2_ll0_temp(keptT))
%             herdTypes(t,1)=3;
%         elseif  r2_sl_temp(t)<= r2_ll_temp(t) & r2_ll0_temp(t)>median(r2_ll0_temp(keptT))
%             herdTypes(t,1)=4;
%         end
%     end
    w_table=readtable([expdir exp '/sound/wordTimeStamp.csv']);
%     w_table.herdTypes=(nan([size(w_table,1) 1]));
%     i=round(w_table.trmin);
%     w_table.herdTypes(find(i>0 & i<=tn))=herdTypes(i(i>0 & i<=tn));
%     
    load([expdir exp '/sound/text_trVectors_binSize' num2str(binSize)   '.mat'],'text_trVectors');
        textvarNames=text_trVectors.Properties.VariableNames;
   
    for ki=1:length(kws);
        kw=kws{ki};
        %% adjsut temporal alignment
        kw_v=text_trVectors.(kw);
        kw_v=[zeros(text_trVectorShift,1)  ; kw_v(1:(end-text_trVectorShift))];
        kw_vs(:,ki)=kw_v;
    end
    
    
    fsize=[30 9];
    figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
    hold on
    ciplot_claire((squeeze(r2_ll(ri,keptT,:)))',keptT,'r',0.3);
    h = findobj(gca,'Type','line');
    set(h,'linewidth',1.5,'linestyle',':')
    hold on
    ciplot_claire((squeeze(r2_sl(ri,keptT,:)))',keptT,'r',0.3);
    hold on
    ciplot_claire((squeeze(r2_ll0(ri,keptT,:)))',keptT,'b',0.3);
    hold on
    plot((kw_v(keptT)/max(kw_v))*0.5,'color','g');
    
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
    scatter(r2_sl_temp(keptT),r2_ll0_temp(keptT),10,'k','filled','markerfacealpha',.1);
       cols=distinguishable_colors(length(textvarNames));
    for ki=1:length(kws);
        kw_v=kw_vs(:,ki);
        kwvi=(kw_v>quantile(kw_v,0.5));
        scatter(r2_sl_temp(kwvi),r2_ll0_temp(kwvi),10,cols(ismember(textvarNames,kws{ki}),:) ,'filled' ,'markerfacealpha',.5);
    end
    hold off;
    legend(cat(2,{''},kws),'location','northwest');
    legend boxoff
    
    
    set(gca,'YScale', 'log','XScale', 'log')
    xlabel('SL coupling'); ylabel('LL coupling');
    title({[upper(exp(1)) exp(2:end) ', ' rnames{ri} ],['shift ' num2str(text_trVectorShift) ' TR']});
    xlim([min(r2_sl_temp(:))-0.01 max(r2_sl_temp(:))+0.01]);
    ylim([min(r2_ll0_temp(:))-0.01 max(r2_ll0_temp(:))+0.01]);
    set(gca,'fontsize',14);
    xlabel({'Speaker-Listener coupling','(R-squared)'});
    ylabel({'Listener-Listener coupling','(R-squared)'});
end




