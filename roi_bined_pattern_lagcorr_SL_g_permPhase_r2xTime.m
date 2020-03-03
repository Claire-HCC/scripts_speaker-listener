clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:-1};
binSize=10;

for ei=1:4;%1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        for si=1:48;
            
            if exist(([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/perm/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL' num2str(si) '.mat']));
                load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/perm/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL' num2str(si) '.mat'],'r2');
                keptT=find(~isnan(r2(1,:)));
                for ri=1:length(rnames);
                    r2Xtime_perm(ri,si)=corr(r2(ri,keptT)',keptT','type','spearman');
                end
            end
        end
        
        load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  ],'rnames','r2');
        keptT=find(~isnan(r2(1,:)));
        for ri=1:length(rnames);
            r2Xtime(ri,1)=corr(r2(ri,keptT)',keptT','type','spearman');
        end
        
        r2Xtime_z=atanh(r2Xtime);
        r2Xtime_perm_z=atanh(r2Xtime_perm);
        
        for ri=1:length(rnames);
            [~,p(ri,1),~,stats] =ttest(r2Xtime_perm_z(ri,:),r2Xtime_z(ri),'tail','left');
            t(ri,1)=-stats.tstat;
        end
        % remove rois showing negative e2Xtime in the permuted data
        [neg_perm_mask]=ttest(r2Xtime_perm_z',0,'tail','left');
        p(neg_perm_mask==1)=NaN;
        t(neg_perm_mask==1)=NaN;
        
        p_fdr=nan(size(p));
        [~,~,p_fdr(~isnan(p))]=fdr(p(~isnan(p)));
        
        table(rnames(p_fdr<.05),r2Xtime(p_fdr<.05),t(p_fdr<.05))
        
        
        save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_r2Xtime' ],'r2Xtime','r2Xtime_perm','p','t','p_fdr','rnames');
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(p_fdr<.05),   r2Xtime(p_fdr<.05)]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_r2Xtime.nii']);
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(p_fdr<.05),   t(p_fdr<.05)]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_r2Xtime_t.nii']);
        
        clear r2Xtime_perm r2Xtime p t p_fdr
    end
end
