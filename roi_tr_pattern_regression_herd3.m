clear all;

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

lags_tested={-10:-4, -20:-4, -30:-4, -10:-1};
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

for ei=3;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_stats'],'r2_sig_fdr','rnames','keptT','r2_byTime_real','r2_byTime_null');
        r2_byTime_s=r2_byTime_real-mean(r2_byTime_null,3);
        SL_sig=r2_sig_fdr;
        
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag0-0_stats' ],'r2_byTime','r2_sig_fdr');
        listenerN=size(r2_byTime,3);
        r2_byTime_l=r2_byTime;
        r2_byTime_lm=mean(r2_byTime,3);
        LL_sig=r2_sig_fdr;
        
        [herd(:,1)]=corr_col(r2_byTime_s(:,keptT)',r2_byTime_lm(:,keptT)');
        
        for perm=1:listenerN;
            r2_byTime_null_temp=r2_byTime_null;
            r2_byTime_null_temp(:,:,perm)=r2_byTime_real;
            r2_byTime_s_temp=r2_byTime_null(:,:,perm)-mean(r2_byTime_null_temp,3);
            
            load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/perm/regression_LL_lag0-0' '_perm' num2str(perm)],'r2_byTime');
            r2_byTime_lm_temp=mean(r2_byTime,3);
            
            herd_null(:,perm)=corr_col(r2_byTime_s_temp(:,keptT)',r2_byTime_lm_temp(:,keptT)');
        end
        herd_z=real(atanh(herd));
        herd_null_z=real(atanh(herd_null));
        
        for ri=1:length(rnames);
            [~,p(ri,1),~,stats]=ttest(herd_null_z(ri,:),herd_z(ri));
            t(ri,1)=-stats.tstat;
        end
        % p=sum(herd_null>herd,2)/listenerN;
        
        % test herding effect within rois showing significant SL and LL
        mask=( SL_sig==1 & LL_sig==1);
        
        sig_fdr_pos(mask,1)=(fdr0(p(mask),0.05)==1 & t(mask)>0);
        sig_fdr_neg(mask,1)=(fdr0(p(mask),0.05)==1 & t(mask)<0);
        
        sig_fwe_pos(mask,1)=(p(mask)<(0.05/sum(mask)) & t(mask)>0);
        sig_fwe_neg(mask,1)=(p(mask)<(0.05/sum(mask)) & t(mask)<0);
        
        save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_regression_SL_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'rnames','lags','herd','herd_null','sig_fdr_pos','sig_fdr_neg','sig_fwe_pos','sig_fwe_neg','p','t');
        
        roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(sig_fdr_pos==1),'UniformOutput',0);
        roi_ids=cell2mat(roi_table.id(cell2mat(roi_table_inds)));
        nii=roiTable2wholeBrainNii_mor([roi_ids, herd(sig_fdr_pos==1)]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_regression' num2str(min(lags)) '-' num2str(max(lags)) '_fdr.nii']);
        
        roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(sig_fdr_neg==1),'UniformOutput',0);
        roi_ids=cell2mat(roi_table.id(cell2mat(roi_table_inds)));
        nii=roiTable2wholeBrainNii_mor([roi_ids, herd(sig_fdr_neg==1)]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_regression' num2str(min(lags)) '-' num2str(max(lags)) '_fdr_neg.nii']);
        
        roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(sig_fwe_pos==1),'UniformOutput',0);
        roi_ids=cell2mat(roi_table.id(cell2mat(roi_table_inds)));
        nii=roiTable2wholeBrainNii_mor([roi_ids, herd(sig_fwe_pos==1)]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_regression' num2str(min(lags)) '-' num2str(max(lags)) '_fwe.nii']);
        
        roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(sig_fwe_neg==1),'UniformOutput',0);
        roi_ids=cell2mat(roi_table.id(cell2mat(roi_table_inds)));
        nii=roiTable2wholeBrainNii_mor([roi_ids, herd(sig_fwe_neg==1)]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_regression' num2str(min(lags)) '-' num2str(max(lags)) '_fwe_neg.nii']);
    end
    
end
