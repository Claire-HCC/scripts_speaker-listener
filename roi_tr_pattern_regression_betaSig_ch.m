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
for ei=1:4;
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
    
    % b_sig_ll
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F','b','rnames');
    b_ll=b(:,2:end,:);
    [~,p]=ttest(b(:,2:end,:),0,'Dim',3);
    b_sig=zeros(size(p));
    p_temp=p(couple_sig_fdr==1 ,:);
    b_sig_ll=zeros(size(p));
    b_sig_ll(couple_sig_fdr==1 ,:)=reshape(fdr0(p_temp(:),0.05),sum(couple_sig_fdr),size(p,2));

  

    % ris=find(couple_sig_fdr==1);
    ris=find(ismember(rnames,{'vPCUN'}));
    for i=1:length(ris);
        ri=ris(i);
        rname=rnames{ri};
        
        figure(fig(i))
        subplot(2,4,ei);
        plot(lags*1.5,squeeze(b_ll(ri,:,:)));
        xlim([-15 15]);
        %   ylim([-0.02 0.12])
             set(gca,'xtick',[-10:2:10],'xtick',[-10:2:10]*1.5);
        grid on
        if sum(b_sig_ll(ri,:))>0;
            text(lags(b_sig_ll(ri,:)==1)*1.5,max(b_ll(ri,b_sig_ll(ri,:)==1,:),[],3),'*','fontsize',20,'HorizontalAlignment','center');
        end
        title(exp);
        
        
        subplot(2,4,ei+4);
        plot(lags*1.5,b_sl(ri,:),'k','linewidth',1.5);
           xlim([-15 15]);
        
        set(gca,'xtick',[-10:2:10],'xtick',[-10:2:10]*1.5);
        %   ylabel('Weightings');
        
        
        if sum(b_sig_sl(ri,:))>0;
            text(lags(b_sig_sl(ri,:)==1)*1.5,b_sl(ri,b_sig_sl(ri,:)==1,1),'*','fontsize',20,'HorizontalAlignment','center');
        end
        %  title({exp,rname,'Speaker-Listener'});
        grid on
       
    end
    
end
