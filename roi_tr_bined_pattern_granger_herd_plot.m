% bined r2 vs std between subjects
clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
rname='vPCUN'
% -4 for merlin, -3 for sherlock
lags=-60:-4;
binSize=30;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

% load([expdir '/scripts_speaker-listener/merlin_' rname '_hmm_findListenersEventInSpeaker.mat'], 'segments_LG');

% for i=(1-min(lags)):length(eventLabels_LG);
%     if length(unique(eventLabels_LG(i+lags)))~=1
%         eventLabels_LG(i)=0;
%     end
% end

for ei=3%:4;%1:4;
    exp=experiments{ei};
    
    eventLabels_LG=h5read([expdir   exp '/fmri/hmm/' rname  '_findListenersEventInSpeaker.hdf5'],'///eventLabels_LG');
    K=length(unique(eventLabels_LG(:,1)));
    cols=distinguishable_colors(K);
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_binsize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2_l','r2_ls','keptT','rnames');
    r2_sl=r2_ls-r2_l;
    keptT_sl=keptT;
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_binSize' num2str(binSize) '_lag0-0' ],'r2','keptT');
    r2_ll=nanmean(r2,4);
    keptT_ll=keptT;
    
    %     load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_lag0-0'],'r2');
    %     r2_ll=nanmean(r2,3);
    %     r2_ll=repmat(r2_ll,1,1,size(r2_sl,3));
    
    % try dividing the story into early and late parts
    keptT=max(min(keptT_ll),min(keptT_sl)):size(r2_sl,2);
    
    eventLabels_LG=eventLabels_LG(keptT);
    
    %     load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'lags','herd','herd_p','herd_sig_fdr_pos','herd_sig_fdr_neg','herd_null');
    %     herdm=mean(herd,2);
    %     herdm_null=nanmean(herd_null,2);
    % rnames(herd_sig_fdr_pos==1)
    
    % ris=find(ismember(rnames,{'HG_L','vPCUN','aCUN'}));
    ris=find(ismember(rnames,rname));
    % [~,ri]=max(herdm);
    
    for rii=1:length(ris);
        ri=ris(rii);
        
        fsize=[20 9];
        figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
        
        r2_sl_temp=r2_sl(:,keptT);
        r2_ll_g=nanmean(r2_ll(:,keptT,:),3);
        
        subplot(1,2,1);
        hold on
        for Ki=1:K;
            a=area(eventLabels_LG==Ki);
            set(a,'facecolor',cols(Ki,:),'linestyle','none','facealpha',0.3);
        end
        p1=plot([r2_sl_temp(ri,:)'],'r','linewidth',2);
        plot([r2_ll_g(ri,:)'],'k','linewidth',2);
        
        % pause(0.01); %  somehow the coloring doesn't work without pause
        % set(p1.Edge, 'ColorBinding','interpolated', 'ColorData', [uint8(c(1:length(keptT),:)*255) uint8(ones(length(keptT),1)*255)].');
        % set(p1.Edge, 'ColorBinding','interpolated', 'ColorData', [uint8(cols(eventLabels_LG,:)*255) uint8(ones(length(keptT),1)*255)].');
        legend({'SL coupling','LL coupling'},'fontsize',13); legend boxoff
        title({exp,rnames{ri}});
        xlim([1 length(keptT)]);
        % ylim([0 1.1*max(max(squeeze(r2_slf(ri,:,:))))]);
        xlabel('Time (TR)'); % ylabel('Coupling');
        set(gca,'fontsize',13);
        
        subplot(1,2,2);
        scatter(r2_sl_temp(ri,:)',r2_ll_g(ri,:)',60,cols(eventLabels_LG,:),'filled');
        [r p]=corr(r2_sl_temp(ri,:)',r2_ll_g(ri,:)');
        xlabel('SL coupling'); ylabel('LL coupling');
        
        % if herd_sig_fdr_pos(ri)==1; star='*'; else  star='';end
        % title(sprintf('Averaged R = %.02f%s',herdm(ri),star));
        
        title(sprintf('Averaged R = %.02f%s',r));
        set(gca,'fontsize',13);
    end
    
    clear herd herd_null herd_p herd_sig_fdr_pos herd_sig_fdr_neg
end
