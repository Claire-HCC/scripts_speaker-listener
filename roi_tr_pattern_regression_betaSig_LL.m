clear all;
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
lags=-40:40;

exp=experiments{2};
load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_exps_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'b_sl','b_ll','keptT','rnames','lags','herd','herd_p','herd_sig');

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');


% ris=find(sum(herd_sig==0,2)==0);
for ei=3:4;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F','b','rnames');
    
    [~,p]=ttest(b(:,2:end,:),0,'Dim',3);
    b_sig=zeros(size(p));
    p_temp=p(herd_sig(:,ei)==1,:);
    b_sig=zeros(size(p));
    b_sig(herd_sig(:,ei)==1,:)=reshape(fdr0(p_temp(:),0.05),sum(herd_sig(:,ei)),size(p,2));
    [bSig_ri bSig_lagi]=find(b_sig);
    figure;
    hist(lags(bSig_lagi));
    
    ris=find(herd_sig(:,1));
    for i=1:length(ris);
        ri=ris(i);
        rname=rnames{ri};
        fsize=[40 20];
        figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
        
        subplot(1,2,1);
        plot(lags,squeeze(b_ll(ri,:,ei,:)));
        grid on
         if sum(b_sig(ri,:))>0;
        text(lags(b_sig(ri,:)==1),b_sl(ri,b_sig(ri,:)==1,ei)+0.001,'*','fontsize',20,'HorizontalAlignment','center');
        end
        title({rname,exp});
        ylabel('Listener to Listener Weightings')
        
        subplot(1,2,2);
        plot(lags,b_sl(ri,:,ei),'k','linewidth',1.5);
        ylabel('Speaker to Listener Weightings');
        
       
        grid on
    end
    xlabel('Speaker precedes-------Time Shift (TR)--------Listeners precede');
end
    
    
    
