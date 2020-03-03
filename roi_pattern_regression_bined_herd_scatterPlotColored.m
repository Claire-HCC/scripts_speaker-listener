
clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
lags=-10:-4;
binSize=10;
rname='vPCUN'
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
exp=experiments{2};
% load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_exps_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'b_sl','b_ll','keptT','rnames','lags','herd','herd_p','herd_sig');

fsize=[35 20];
fig1=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
fig2=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);

for ei=1%1:4;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/LLselfother/regression_LL_binSize' num2str(binSize) '_lag0-0' ],'r2','rnames');
    r2_ll=nanmean(r2,3);
    
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
    r2_sl=r2;
    
    keptT=[min(find(nansum(r2_sl,1)>0)):max(find(nansum(r2_sl)>0))];
    
    ri=find(ismember(rnames,rname));
    
    r2_ll_temp=r2_ll(ri,keptT);
    r2_sl_temp=r2_sl(ri,keptT);
    
    cols=zeros(length(keptT),3);
    colx=r2_sl_temp>median(r2_sl_temp);%((max(r2_sl_temp)-min(r2_sl_temp))/2+min(r2_sl_temp));
    coly=r2_ll_temp>median(r2_ll_temp);%((max(r2_ll_temp)-min(r2_ll_temp))/2+min(r2_ll_temp));
    
    
    set(0, 'currentfigure', fig1);
    for ti=1:length(keptT);
        %    coli=1+ceil(atan2d(coupling_sl_temp(ti)>median(coupling_sl_temp),r2_ll_temp(ti)-median(r2_ll_temp)));
        if  colx(ti)==1 & coly(ti)==1;
            cols(ti,:)=[1 0.7 0];
        elseif colx(ti)==1 & coly(ti)==0;
            cols(ti,:)=[0 1 0];
        elseif colx(ti)==0 & coly(ti)==0;
            cols(ti,:)=[0 0 1];
        elseif colx(ti)==0 & coly(ti)>0;
            cols(ti,:)=[0.7 0.1 1];
        end
    end
    
    set(0, 'currentfigure', fig1);
    hold on
    for ti=1:length(keptT);
        %  x=(keptT(ti)-binSize/2+1):(keptT(ti)+binSize/2);
        %   area(x,ones(1,binSize),'Facecolor',cols(ti,:),'FaceAlpha',0.3,'linestyle','none');
        
        x=keptT(ti);
        line([x x],[0 1],'color',[cols(ti,:) 0.3])
    end
    
    p1=plot(keptT,[r2_sl_temp],'r','linewidth',2);
    plot(keptT,[r2_ll_temp],'b','linewidth',2);
    
    hold off
    title({[exp ', ' rnames{ri}],['binSize' num2str(binSize) ', Lags:' num2str(min(lags)) '~' num2str(max(lags)) ]});
    xlim([min(keptT)-binSize/2 max(keptT)+binSize/2]);
    ylim([0 max(r2_sl_temp(:))+0.05]);
    
    % ylim([0 1.1*max(max(squeeze(coupling_slf(ri,:,:))))]);
    xlabel('Time (TR)'); % ylabel('Coupling');
    set(gca,'fontsize',13);
    
    
    set(0, 'currentfigure', fig2);
    hold on
    %        subplot(length(binSize_tested),length(lags_tested),(bi-1)*3+lagi);
    for ti=1:length(keptT);
        coli=1+ceil(atan2d(r2_sl_temp(ti)-median(r2_sl_temp),r2_ll_temp(ti)-median(r2_ll_temp)));
        if coli<0; coli=360+coli;end
        scatter(r2_sl_temp(ti),r2_ll_temp(ti),60,cols(ti,:),'filled');
    end
    
    % scatter(r2_sl_temp,r2_ll_temp,60,'k','filled');hold on
    % scatter(r2_sl_temp(p_temp==0),r2_ll_temp(ri,p_temp==0)',60,'b','filled');
    [r ]=corr(r2_sl_temp,r2_ll_temp);
    xlabel('SL coupling'); ylabel('LL coupling');
    hold off
    % if herd_sig_fdr_pos(ri)==1; star='*'; else  star='';end
    % title(sprintf('Averaged R = %.02f%s',herdm(ri),star));
    title({exp,rname})
    %    title({sprintf('Averaged R = %.2f, p = %.3f',r,p_herd(ri)),['binSize' num2str(binSize) ', Lags:' num2str(min(lags)) '~' num2str(max(lags)) ]});
    set(gca,'fontsize',13);
end



