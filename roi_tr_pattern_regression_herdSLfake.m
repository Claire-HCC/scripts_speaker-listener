% bined r2 vs std between subjects
clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

% -4 for merlin, -3 for sherlock
lags=-6:-3;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

for ei=1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_SL_minus1listener_lag' num2str(min(lags)) '-' num2str(max(lags))],'r2','rnames','keptT');
    r2_sl=r2;
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_SLfake_lag' num2str(min(lags)) '-' num2str(max(lags))],'r2');
    r2_slf=r2;
    
    % load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_minus1listener_lag0-0'],'r2');
    %     r2_ll=nanmean(r2,4);
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_lag0-0'],'r2');
    r2_ll=nanmean(r2,3);
    r2_ll=repmat(r2_ll,1,1,size(r2_sl,3));
    
    % try dividing the story into early and late parts
    keptT=min(keptT):size(r2_sl,2);
    
    for sp=1:size(r2,3);
        for ri=1:size(r2,1)'
            herd(ri,sp)=corr(r2_sl(ri,keptT,sp)',r2_ll(ri,keptT,sp)','type','spearman');
            herd_null(ri,sp)=corr(r2_slf(ri,keptT,sp)',r2_ll(ri,keptT,sp)','type','spearman');
        end
    end
    
    herdz=real(atanh(herd));
    herdz_m=mean(herdz,2);
    herdz_null=real(atanh(herd_null));
    d=(herdz-herdz_null)';
    d=mean(d);
    [~,herd_p]=ttest([herdz-herdz_null]',0);
    
     herd_sig_fdr_pos=zeros(size(herd_p));
    herd_sig_fdr_pos(d>0)=fdr0(herd_p(d>0),0.025);
         herd_sig_fdr_neg=zeros(size(herd_p));
      herd_sig_fdr_neg(d<0)=fdr0(herd_p(d<0),0.025);

    % test herd effect within regions showing significant coupling
    % load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag-10-10_stats'],'couple_sig_fdr');
    % herd_sig_fdr(couple_sig_fdr==1)=fdr0(herd_p(couple_sig_fdr==1),0.05);
    
    herdm=mean(herd,2);
    herdm_null=nanmean(herd_null,2);
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames);
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, herdm]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '.nii']);
    
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(herd_sig_fdr_pos==1));
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, herdm(herd_sig_fdr_pos==1)]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr_pos.nii']);
    
    roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(herd_sig_fdr_neg==1));
    roi_ids=cell2mat(roi_table.id(roi_table_inds));
    nii=roiTable2wholeBrainNii_mor([roi_ids, herdm(herd_sig_fdr_neg==1)]);
    save_nii(nii,[expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '_fdr_neg.nii']);
    
    rnames(herd_sig_fdr_pos==1)
    save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'lags','herd','herd_p','herd_sig_fdr_pos','herd_sig_fdr_neg','herd_null');
    
    
    c=hot(length(keptT)+100);
    ris=find(ismember(rnames,{'HG_L','vPCUN','aCUN'}));
    % [~,ri]=max(herdm);
    for rii=1:length(ris);
        ri=ris(rii);
        figure;
        for sp=1:size(r2,3);
            r2_sl_g=nanmean(r2_sl(:,keptT,sp),3);
            r2_ll_g=nanmean(r2_ll(:,keptT,sp),3);
            r2_slf_g=nanmean(r2_slf(:,keptT,sp),3);
            ylim([min([r2_sl_g(ri,:) r2_slf_g(ri,:)]) max([r2_sl_g(ri,:) r2_slf_g(ri,:)])]);
            
            subplot(2,2,1);
            scatter(r2_sl_g(ri,:)',r2_ll_g(ri,:)',50,c(1:length(keptT),:),'filled');
            xlabel('SL coupling'); ylabel('LL coupling');
            xlim([min([r2_sl_g(ri,:) r2_slf_g(ri,:)]) max([r2_sl_g(ri,:) r2_slf_g(ri,:)])]);
            if herd_sig_fdr_pos(ri)==1; star='*'; else  star='';end
            title(sprintf('herd R=%.02f%s',herdm(ri),star));
            hold on
            
            subplot(2,2,2);
            p1=plot([r2_sl_g(ri,:)'],'r');
            
            hold on
            plot([r2_ll_g(ri,:)'],'b','linewidth',2);
            pause(0.01); %  somehow the coloring doesn't work without pause
            set(p1.Edge, 'ColorBinding','interpolated', 'ColorData', [uint8(c(1:length(keptT),:)*255) uint8(ones(length(keptT),1))].');
            
            legend('SL coupling','LL coupling'); legend boxoff
            title(rnames{ri})
            xlim([1 length(keptT)]);
            ylim([0 max([r2_sl_g(ri,:) r2_slf_g(ri,:) ])+0.05]);
            
            subplot(2,2,3);
            scatter(r2_slf_g(ri,:)',r2_ll_g(ri,:)',50,c(1:length(keptT),:),'filled');
            xlabel('pseudo-SL coupling'); ylabel('LL coupling');
            xlim([min([r2_sl_g(ri,:) r2_slf_g(ri,:)]) max([r2_sl_g(ri,:) r2_slf_g(ri,:)])]);
            title(sprintf('pseudo-herd R=%.02f',herdm_null(ri)));
            hold on
            
            subplot(2,2,4);
            p2=plot([r2_slf_g(ri,:)'],'r');
            pause(0.01)
            set(p2.Edge, 'ColorBinding','interpolated', 'ColorData', [uint8(c(1:length(keptT),:)*255) uint8(ones(length(keptT),1))].');
            hold on
            plot([r2_ll_g(ri,:)'],'b','linewidth',2);
            legend('pseudo-SL coupling','LL coupling'); legend boxoff
            title(rnames{ri});
             xlim([1 length(keptT)]);
             ylim([0 max([r2_sl_g(ri,:) r2_slf_g(ri,:) ])+0.05]);
            
            hold on
            
            % clear p1 p2
        end
    end
    
    
    clear herd herd_null herd_p herd_sig_fdr_pos herd_sig_fdr_neg
end
