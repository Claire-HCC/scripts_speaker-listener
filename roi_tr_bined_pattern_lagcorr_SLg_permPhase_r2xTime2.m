clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:-4};
binSize=10;

for ei=1:4;%1:4;%1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        for si=1:48;
            if exist(([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/perm/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL' num2str(si) '.mat']));
                load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/perm/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_permSL' num2str(si) '.mat'],'r2');
                keptT=find(~isnan(r2(1,:)));
                r2_perm(:,:,si)=r2;
                for ri=1:length(rnames);
                    
                    r2Xtime_perm(ri,si)=corr(r2(ri,keptT)',keptT','type','spearman');
                end
                
                load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/perm/regression_SL_binSize' num2str(binSize) '_lag-3-3_permSL' num2str(si) '.mat'],'r2');
                r2_perm0(:,:,si)=r2;
                keptT=find(~isnan(r2(1,:)));
                for ri=1:length(rnames);
                    r2Xtime_perm0(ri,si)=corr(r2(ri,keptT)',keptT','type','spearman');
                end
            else
                r2Xtime_perm0(:,si)=NaN;
                r2Xtime_perm(:,si)=NaN;
                r2_perm(:,:,si)=NaN;
                r2_perm0(:,:,si)=NaN;
            end
        end
        
        
        
        load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/regression_SL_binSize' num2str(binSize) '_lag-3-3'  ],'rnames','r2');
        keptT=find(~isnan(r2(1,:)));
        r20=r2;
        for ri=1:length(rnames);
            r2Xtime0(ri,1)=corr(r2(ri,keptT)',keptT','type','spearman');
        end
        
        load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  ],'rnames','r2');
        
        keptT=find(~isnan(r2(1,:)));
        for ri=1:length(rnames);
            r2Xtime(ri,1)=corr(r2(ri,keptT)',keptT','type','spearman');
        end
        
        listenerN=size(r2_perm,3)
        d_real=repmat(atanh(r2Xtime),1,listenerN)-atanh(r2Xtime_perm);
         d_null=repmat(atanh(r2Xtime0),1,listenerN)-atanh(r2Xtime_perm0);
        
        [~,p,~,stats]= ttest(d_real',d_null','tail','right');
        p=p';
        t=stats.tstat';
%         clear p t
%         for ri=1:length(rnames);
%             [~,p(ri,1),~,stats] =ttest(d_null(ri,:),d_real(ri),'tail','left');
%             t(ri,1)=-stats.tstat;
%         end
        
           [neg_perm_mask]=ttest(d_null',0,'tail','left');
        p(neg_perm_mask==1)=NaN;
        t(neg_perm_mask==1)=NaN;
        
        
        %         % remove rois showing negative e2Xtime in the permuted data
        %         [neg_perm_mask]=ttest(r2Xtime_perm_z',0,'tail','left');
        %         p(neg_perm_mask==1)=NaN;
        %         t(neg_perm_mask==1)=NaN;
        
        p_fdr=nan(size(p));
        [~,~,p_fdr(~isnan(p))]=fdr(p(~isnan(p)));
        
        table(rnames(p_fdr<.05),r2Xtime(p_fdr<.05),t(p_fdr<.05))
        
        
%         save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_r2Xtime' ],'r2Xtime','r2Xtime_perm','p','t','p_fdr','rnames','r2','r20','r2_perm','r2_perm0');
%         %
%         nii=roiTable2wholeBrainNii_mor([roi_ids(p_fdr<.05),   r2Xtime(p_fdr<.05)]);
%         save_nii(nii,[expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_r2Xtime.nii']);
%         
        %         nii=roiTable2wholeBrainNii_mor([roi_ids(p_fdr<.05),   t(p_fdr<.05)]);
        %         save_nii(nii,[expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_r2Xtime_t.nii']);
        %
        clear r2Xtime_perm r2Xtime p t p_fdr clear r2 r2_perm r20 r2_perm0
    end
end
