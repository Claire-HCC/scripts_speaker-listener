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

lags_tested={-4:4,-15:15,-10:10,  -40:40, -10:-1 , -40:-1, 1:10,  1:40 };

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase.mat' ],'b','r2','rnames');
        r2_perm=r2;
        b_perm=b;
        
        load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat' ],'b','r2','rnames');
        
        p=mean(r2_perm>r2,2);
        p(isnan(r2))=NaN;
        
        pfdr=nan(size(p));
        [~,~,pfdr(~isnan(p))]=fdr(p(~isnan(p)));
        sig_fdr=pfdr<.05;
        
        p_fwemax=sum(repmat(max(r2_perm),length(r2),1)>r2,2)/(size(r2_perm,2));;
        p_fwemax(isnan(r2))=NaN;
        
        %     subplot(4,1,ei);
        %     boxplot(r2_perm');
        %     hold on
        %     scatter(1:61,r2,40,'r','filled')
        %     set(gca,'xtick',1:61,'xticklabels',rnames);
        %     xtickangle(45)
        %     xlim([0 61]);
        %     title(exp)
        %     grid on
        %     text(find(pfdr<.05),r2(pfdr<.05),'*','fontsize',14);
        
        save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'r2','sig_fdr','p','pfdr','p_fwemax','r2_perm','b_perm');
        %
        rnames(pfdr<.05)
        %     nii=roiTable2wholeBrainNii_mor([roi_ids(pfdr<.05),  r2(pfdr<.05)]);
        %     save_nii(nii,[expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_fdr.nii']);
        %
      %  rnames(p_fwemax<.05)
        %     nii=roiTable2wholeBrainNii_mor([roi_ids( p_fwemax<.05),  r2( p_fwemax<.05)]);
        %     save_nii(nii,[expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_g/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_ p_fwemax.nii']);
    end
end



