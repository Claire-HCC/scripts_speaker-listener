% bined r2 vs std between subjects
clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
K=28;

% -4 for merlin, -3 for sherlock
lags=-7:-4;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

load([expdir '/scripts_speaker-listener/merlin_dPCC_hmm_findListenersEventInSpeaker.mat'], 'segments_LG');
cols=distinguishable_colors(K);
segments_LG=segments_LG(K,:);

for i=(1-min(lags)):length(segments_LG);
    if length(unique(segments_LG(i+lags)))~=1
        segments_LG(i)=0;
    end
end


for ei=3%:4;%1:4;
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
    
    segments_LG=segments_LG(keptT);
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'lags','herd','herd_p','herd_sig_fdr_pos','herd_sig_fdr_neg','herd_null');
    herdm=mean(herd,2);
    herdm_null=nanmean(herd_null,2);
    rnames(herd_sig_fdr_pos==1)
    
    % ris=find(ismember(rnames,{'HG_L','vPCUN','aCUN'}));
    ris=find(ismember(rnames,{'aCUN'}));
    % [~,ri]=max(herdm);
    
    for rii=1:length(ris);
        ri=ris(rii);
        
        fsize=[20 9];
        figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
        
        r2_sl_g=nanmean(r2_sl(:,keptT,:),3);
        r2_ll_g=nanmean(r2_ll(:,keptT,:),3);
        r2_slf_g=nanmean(r2_slf(:,keptT,:),3);
        
        subplot(1,2,1);
        p1=plot([r2_sl_g(ri,:)'],'r','linewidth',2);
        hold on
        plot([r2_ll_g(ri,:)'],'k','linewidth',2);
        pause(0.01); %  somehow the coloring doesn't work without pause
        % set(p1.Edge, 'ColorBinding','interpolated', 'ColorData', [uint8(c(1:length(keptT),:)*255) uint8(ones(length(keptT),1)*255)].');
        set(p1.Edge, 'ColorBinding','interpolated', 'ColorData', [uint8(cols(segments_LG,:)*255) uint8(ones(length(keptT),1)*255)].');
        legend({'SL coupling','LL coupling'},'fontsize',13); legend boxoff
        title({exp,rnames{ri}});
        xlim([1 length(keptT)]);
        % ylim([0 1.1*max(max(squeeze(r2_slf(ri,:,:))))]);
        xlabel('Time (TR)'); % ylabel('Coupling');
        set(gca,'fontsize',13);
        
        subplot(1,2,2);
        scatter(r2_sl_g(ri,:)',r2_ll_g(ri,:)',60,cols(segments_LG,:),'filled');
        
        xlabel('SL coupling'); ylabel('LL coupling');
        
        if herd_sig_fdr_pos(ri)==1; star='*'; else  star='';end
        title(sprintf('Averaged R = %.02f%s',herdm(ri),star));
        set(gca,'fontsize',13);
    end
    
    clear herd herd_null herd_p herd_sig_fdr_pos herd_sig_fdr_neg
end
