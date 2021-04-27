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
lags_tested={-10:-4. 4:10, -10:10, -10:-1, 0, 1:10};

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_each/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase.mat' ],'b','r2','rnames');
        r2_perm=r2;
        b_perm=b;
        load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat' ],'b','r2','rnames');
        r2=r2;
        b=b;
        
        p=mean(nanmean(r2_perm,3)>nanmean(r2,3),4);
        p(isnan(r2(:,1,1)))=NaN;
        pfdr=nan(size(p));
        [~,~,pfdr(~isnan(p(:)))]=fdr(p(~isnan(p)));
        pfdr=reshape(pfdr,size(p));
        
        pfwe=p*sum(~isnan(p(:)));
        
        save([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'r2','p','pfdr','pfwe','b');
        
        rnames(pfdr<.05)
        nii=roiTable2wholeBrainNii_mor([roi_ids(pfdr<.05),  r2(pfdr<.05)]);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_each/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_pfdr.nii']);
        
        %      rnames(pfwe<.05)
        nii=roiTable2wholeBrainNii_mor([roi_ids( pfwe<.05),  r2( pfwe<.05)]);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_each/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_ pfwe.nii']);
    end
end


