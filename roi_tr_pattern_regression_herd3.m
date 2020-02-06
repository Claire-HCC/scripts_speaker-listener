clear all;
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

lags_tested={-10:-4, -20:-4, -30:-4, -10:-1};
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
roi_ids=cell2mat(roi_table.id);
rnames=table2array(roi_table(:,3));

for ei=1:4%:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/SLeach/regression_SLeach_lag-10-10_classification' ],'sig_fdr_withinSigSL');
        sig_betaClass=sig_fdr_withinSigSL;
        
        load(sprintf('%s/%s/fmri/pattern_regression/%s/roi/%s/LLselfother/regression_LL_lag0-0',expdir,exp,timeUnit,froidir),'rnames','r2_byTime');
        r2_byTime_llm=nanmean(r2_byTime,3);
        
        load(sprintf('%s/%s/fmri/pattern_regression/%s/roi/%s/SLg/regression_SL_lag%d-%d',expdir,exp,timeUnit,froidir,min(lags),max(lags)),'r2_byTime','rnames','keptT');%,'r2_sig_fdr');
        r2_byTime_sl=r2_byTime;
        
        [herd(:,1)]=corr_col(r2_byTime_sl(:,keptT)',r2_byTime_llm(:,keptT)');
        
        for perm=1:48;
            if exist(sprintf('%s/%s/fmri/pattern_regression/%s/roi/%s/LLselfother/perm/regression_LL_lag0-0_permSL%03d.mat',expdir,exp,timeUnit,froidir,perm));
                load(sprintf('%s/%s/fmri/pattern_regression/%s/roi/%s/SLg/perm/regression_SL_lag%d-%d_permSL%03d',expdir,exp,timeUnit,froidir,min(lags),max(lags),perm),'rnames','keptT','r2_byTime');
                r2_byTime_sl=r2_byTime;
                
                load(sprintf('%s/%s/fmri/pattern_regression/%s/roi/%s/LLselfother/perm/regression_LL_lag0-0_permSL%03d',expdir,exp,timeUnit,froidir,perm),'r2_byTime','rnames');
                r2_byTime_llm=nanmean(r2_byTime,3);
                
                herd_null(:,perm)=corr_col(r2_byTime_sl(:,keptT)',r2_byTime_llm(:,keptT)');
            end
        end
        
        herd_z=real(atanh(herd));
        herd_null_z=real(atanh(herd_null));
        
        for ri=1:length(rnames);
            rname=rnames(ri);
            [~,p(ri,1),~,stats]=ttest(herd_null_z(ri,:),herd_z(ri),'tail','left');
            t(ri,1)=-stats.tstat;
        end
        
        p_adj_withinBetaClass=nan(size(p));
        [~,~,p_adj_withinBetaClass(sig_betaClass)]=fdr(p(sig_betaClass));
        exp
        min(lags)
        rnames(p_adj_withinBetaClass<.05)
        
        
        save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd/herd_regression_SL_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'rnames','lags','herd','herd_null','p_adj_withinBetaClass','t');
        
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(p_adj_withinBetaClass<.05), t(p_adj_withinBetaClass<.05)]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd/herd_regression' num2str(min(lags)) '-' num2str(max(lags)) '_fdr_withinBetaClass.nii']);
        
        
        figure
        boxplot(herd_null_z');
        % scatter(repmat(1:61,1,48),herd_null_z(:),20,'k','filled','MarkerFaceAlpha',.1);
        hold on
        scatter(1:61,herd_z,40,'r','filled')
        set(gca,'xtick',1:61,'xticklabels',rnames);
        xtickangle(45)
        xlim([0 61]);
        ylabel(exp)
        grid on
        text(find(p_adj_withinBetaClass<.05),herd_z(p_adj_withinBetaClass<.05),'*','fontsize',14);
        title(exp)
        
        clear herd_z herd_null_z herd herd_null
    end
    
end

