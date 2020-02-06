clear all
close all

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
roi_ids=cell2mat(roi_table.id);
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;

lags=-10:10;

figure
for ei=1:4;
    exp=experiments{ei};
    
    if exist(sprintf('%s/%s/fmri/temporal_regression/%s/roi/%s/SLg/perm/regression_SL_lag%d-%d_permPhase.mat',expdir,exp,timeUnit,froidir,min(lags),max(lags)));
        load(sprintf('%s/%s/fmri/temporal_regression/%s/roi/%s/SLg/perm/regression_SL_lag%d-%d_permPhase.mat',expdir,exp,timeUnit,froidir,min(lags),max(lags)),'b_perm','r2_perm','rnames');
    else
        b_perm=[];        r2_perm=[];
        for perm=1:1000;
            load(sprintf('%s/%s/fmri/temporal_regression/%s/roi/%s/SLg/perm/regression_SL_lag%d-%d_permPhase%04d.mat',expdir,exp,timeUnit,froidir,min(lags),max(lags),perm),'b','r2','rnames');
            b_perm(:,:,perm)=b;
            r2_perm(:,perm)=r2;
        end
        save(sprintf('%s/%s/fmri/temporal_regression/%s/roi/%s/SLg/perm/regression_SL_lag%d-%d_permPhase.mat',expdir,exp,timeUnit,froidir,min(lags),max(lags)),'b_perm','r2_perm','rnames');
    end
    
    load([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat' ],'b','r2','rnames');
    
     p=sum(r2_perm>r2(:,1),2)/size(r2_perm,2);
    p_fdr=nan(size(p));
    ris=find(~isnan(r2(:,1)));
    [~,~,p_fdr(ris)]=fdr(p(ris));
    sig_fdr=p_fdr<.05;
    
    p_fwemax=sum(repmat(max(r2_perm),length(r2),1)>r2,2)/(size(r2_perm,2));;
    p_fwemax(isnan(r2))=NaN;
    
    %     subplot(4,1,ei);
    %     boxplot(r2_perm');
    %     hold on
    %     scatter(1:61,r2,40,'r','filled')
    %     set(gca,'xtick',1:61,'xticklabels',rnames);
    %     xtickangle(45)
    %     xlim([0 62]);
    %     title(exp)
    %     grid on
    %     text(find(p_adj<.05),r2(p_adj<.05),'*','fontsize',14);
    %     ylim([0 0.4]);
    
    save([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'r2','sig_fdr','p','p_fdr','p_fwemax','r2_perm','b_perm','b');
    
    rnames(p_fdr<.05)
    nii=roiTable2wholeBrainNii_mor([roi_ids(p_fdr<.05),  r2(p_fdr<.05)]);
    save_nii(nii,[expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_fdr.nii']);
    
    rnames(p_fwemax<.05)
    nii=roiTable2wholeBrainNii_mor([roi_ids( p_fwemax<.05),  r2( p_fwemax<.05)]);
    save_nii(nii,[expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_ p_fwemax.nii']);
end


