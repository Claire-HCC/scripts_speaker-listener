clear all;
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
lags=-10:10;

% ris=find(sum(herd_sig==0,2)==0);
fsize=[30 18];
fig(1)=figure('unit','centimeter','position',[1 1 fsize],'paperposition',[1 1 fsize],'papersize',fsize);
fig(2)=figure('unit','centimeter','position',[1 1 fsize],'paperposition',[1 1 fsize],'papersize',fsize);

fsize=[30 18];
figscatter=figure('unit','centimeter','position',[1 1 fsize],'paperposition',[1 1 fsize],'papersize',fsize);
for ei=3;%1:4;
    exp=experiments{ei};
    
    % b_sig_sl
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F','b','rnames');
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_statsBeta' ],'b_sig');
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_stats' ],'couple_sig_fdr');
    b_sig_sl=b_sig;
    b_sl=b(:,2:end);
    
    % find significant betas withing regions showing significant sl
    % coupling
    [bSig_sl_ri bSig_sl_lagi]=find(b_sig_sl);
    
    figure(figscatter);
    subplot(2,4,ei+4);
    % focus on DMN network
    histogram(lags(bSig_sl_lagi),round(length(lags)/2),'FaceColor',[0.5 0.5 0.5],'binMethod','integers');
    set(gca,'xtick',lags);grid on
    title('Speaker-Listener');
    if ei==4;
        xlabel('Speaker precedes---Lags (TR)---Listeners precede');
    end
    ylabel('Frequency');
    xlim([min(lags)-1 max(lags)+1]);
    ylim([0 35]);
    
    % b_sig_ll
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F','b','rnames');
    b_ll=b(:,2:end,:);
    [~,p]=ttest(b(:,2:end,:),0,'Dim',3);
    b_sig=zeros(size(p));
    p_temp=p(couple_sig_fdr==1 ,:);
    b_sig_ll=zeros(size(p));
    b_sig_ll(couple_sig_fdr==1 ,:)=reshape(fdr0(p_temp(:),0.05),sum(couple_sig_fdr),size(p,2));
    [bSig_ll_ri bSig_ll_lagi]=find(b_sig_ll);
    
    subplot(2,4,ei);
    % focus on DMN network
    histogram(lags(bSig_ll_lagi),round(length(lags)/2),'FaceColor',[0.5 0.5 0.5],'binMethod','integers');
    set(gca,'xtick',lags);grid on
    ylabel('Frequency');
    title({exp,'Listener-Listener'});
    xlabel('Lags (TR)');
    xlim([min(lags)-1 max(lags)+1]);
    %      ylim([0 35]);
    
    
    % ris=find(couple_sig_fdr==1);
    ris=find(ismember(rnames,{'STC_R'}));
    for i=1:length(ris);
        ri=ris(i);
        rname=rnames{ri};
        
        figure(fig(i))
        subplot(2,4,ei);
        plot(lags,squeeze(b_ll(ri,:,:)));
        %   ylim([-0.02 0.12])
        set(gca,'xtick',lags);
        grid on
        if sum(b_sig_ll(ri,:))>0;
            text(lags(b_sig_ll(ri,:)==1),max(b_ll(ri,b_sig_ll(ri,:)==1,:),[],3),'*','fontsize',20,'HorizontalAlignment','center');
        end
        title({exp,rname,'Listener-Listener'});
        ylabel('Weightings');
        xlabel('Lags (TR)');
        
        subplot(2,4,ei+4);
        plot(lags,b_sl(ri,:),'k','linewidth',1.5);
        %     ylim([-0.02 0.02]);
        set(gca,'xtick',lags);
        ylabel('Weightings');
        
        
        if sum(b_sig_sl(ri,:))>0;
            text(lags(b_sig_sl(ri,:)==1),b_sl(ri,b_sig_sl(ri,:)==1,1),'*','fontsize',20,'HorizontalAlignment','center');
        end
        title({exp,rname,'Speaker-Listener'});
        grid on
        if ei==4;
            xlabel('Speaker precedes---Lags (TR)---Listeners precede');
        end
    end
    
end

%
%
%     nnames=unique(roi_table.network);
%     for ni=1:6;
%         nname=nnames{ni};
%
%         fsize=[20 8];
%         figure('unit','centimeter','position',[1 1 fsize],'paperposition',[1 1 fsize],'papersize',fsize);
%
%         % b_sig_sl
%         load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F','b','rnames');
%         load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],'F_perm','b_perm','rnames','couple_sig_fdr');
%          b_sl=b(:,2:end);
%         p=sum(repmat(b(:,2:end),1,1,1000)<b_perm(:,2:end,:),3)/1000;
%         b_sig=zeros(size(p));
%         p_temp=p(couple_sig_fdr==1,:);
%         b_sig_sl=zeros(size(p));
%         b_sig_sl(couple_sig_fdr==1,:)=reshape(fdr0(p_temp(:),0.05),sum(couple_sig_fdr),size(p,2));
%         [bSig_sl_ri bSig_sl_lagi]=find(b_sig_sl);
%         subplot(1,2,2);
%         % focus on DMN network
%         roi_n=ismember(bSig_sl_ri,strmatch(nname,roi_table.network));
%         histogram(lags(bSig_sl_lagi(roi_n)),round(length(lags)/2),'FaceColor',[0.5 0.5 0.5]);
%         title('Speaker-Listener');
%         xlabel('Speaker precedes-------Time Shift (TR)--------Listeners precede');
%         ylabel('Weightings');
%         xlim([min(lags)-1 max(lags)+1]);
%
%         % b_sig_ll
%         load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F','b','rnames');
%         b_ll=b(:,2:end,:);
%         [~,p]=ttest(b(:,2:end,:),0,'Dim',3);
%         b_sig=zeros(size(p));
%         p_temp=p(couple_sig_fdr==1,:);
%         b_sig_ll=zeros(size(p));
%         b_sig_ll(couple_sig_fdr==1,:)=reshape(fdr0(p_temp(:),0.05),sum(couple_sig_fdr),size(p,2));
%         [bSig_ll_ri bSig_ll_lagi]=find(b_sig_ll);
%         subplot(1,2,1);
%         % focus on DMN network
%         roi_n=ismember(bSig_ll_ri,strmatch(nname,roi_table.network));
%         histogram(lags(bSig_ll_lagi(roi_n)),round(length(lags)/2),'FaceColor',[0.5 0.5 0.5]);
%         ylabel('Weightings');
%         title({exp,nname,'Listener-Listener'});
%         xlabel('Time Shift (TR)');
%        xlim([min(lags)-1 max(lags)+1]);
%     end
%
%
