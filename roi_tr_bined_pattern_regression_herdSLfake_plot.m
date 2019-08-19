% bined r2 vs std between subjects
clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

% -4 for merlin, -3 for sherlock
lags=-7:-4;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

for ei=1%:4;%1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_SL_minus1listener_lag' num2str(min(lags)) '-' num2str(max(lags))],'r2','rnames','keptT');
    r2_sl=r2;
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_SLfake_lag' num2str(min(lags)) '-' num2str(max(lags))],'r2');
    r2_slf=r2;
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_minus1listener_lag0-0'],'r2');
        r2_ll=nanmean(r2,4);
%     load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_lag0-0'],'r2');
%     r2_ll=nanmean(r2,3);
%     r2_ll=repmat(r2_ll,1,1,size(r2_sl,3));
    
    % try dividing the story into early and late parts
    keptT=min(keptT):size(r2_sl,2);
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'lags','herd','herd_p','herd_sig_fdr_pos','herd_sig_fdr_neg','herd_null');
    herdm=mean(herd,2);
    herdm_null=nanmean(herd_null,2);
    rnames(herd_sig_fdr_pos==1)
    
    c=hot(length(keptT)+100);
    % ris=find(ismember(rnames,{'HG_L','vPCUN','aCUN'}));
    ris=find(ismember(rnames,{'HG_R','LOC_L','DLPFC_L'}));
    % [~,ri]=max(herdm);
    for rii=1:length(ris);
        ri=ris(rii);
        
        fsize=[34 30];
        figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
        for sp=1:size(r2,3);
            r2_sl_g=nanmean(r2_sl(:,keptT,sp),3);
            r2_ll_g=nanmean(r2_ll(:,keptT,sp),3);
            r2_slf_g=nanmean(r2_slf(:,keptT,sp),3);
            
            
            subplot(2,2,1);
            scatter(r2_sl_g(ri,:)',r2_ll_g(ri,:)',20,c(1:length(keptT),:),'filled');
            xlabel('SL coupling'); ylabel('LL coupling');
      
            if herd_sig_fdr_pos(ri)==1; star='*'; else  star='';end
            title(sprintf('herd R=%.02f%s',herdm(ri),star));
            hold on
            
            subplot(2,2,2);
            p1=plot([r2_sl_g(ri,:)'],'r','linewidth',0.5);
            hold on
            plot([r2_ll_g(ri,:)'],'b','linewidth',2);
            pause(0.01); %  somehow the coloring doesn't work without pause
            set(p1.Edge, 'ColorBinding','interpolated', 'ColorData', [uint8(c(1:length(keptT),:)*255) uint8(ones(length(keptT),1)*255)].');
            legend({'SL coupling','LL coupling'},'fontsize',18); legend boxoff
            title({exp,rnames{ri}});
            xlim([1 length(keptT)]);
            ylim([0 1.1*max(max(squeeze(r2_slf(ri,:,:))))]);
            xlabel('Time (TR)'); ylabel('Coupling');
            
            subplot(2,2,3);
            scatter(r2_slf_g(ri,:)',r2_ll_g(ri,:)',20,c(1:length(keptT),:),'filled');
            xlabel('pseudo-SL coupling'); ylabel('LL coupling');
            title(sprintf('pseudo-herd R=%.02f',herdm_null(ri)));
            hold on
            
            subplot(2,2,4);
            p2=plot([r2_slf_g(ri,:)'],'r','linewidth',0.5);
            pause(0.01)
            set(p2.Edge, 'ColorBinding','interpolated', 'ColorData', [uint8(c(1:length(keptT),:)*255) uint8(ones(length(keptT),1)*255)].');
            hold on
            plot([r2_ll_g(ri,:)'],'b','linewidth',2);
            
            legend({'SL coupling','LL coupling'},'fontsize',18); legend boxoff
            title({exp,rnames{ri}});
            xlim([1 length(keptT)]);
             ylim([0 1.1*max(max(squeeze(r2_slf(ri,:,:))))]);
            xlabel('Time (TR)'); ylabel('Coupling');
            hold on
            
            % clear p1 p2
        end
    end
    
    
    clear herd herd_null herd_p herd_sig_fdr_pos herd_sig_fdr_neg
end
