clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

binSize_tested=[10 15 30]; % tr;
lags_tested={ -10:-4, -20:-4, -30:-4. -10:-1};
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
rnames_table=table2array(roi_table(:,3));
roi_ids=cell2mat(table2array(roi_table(:,1)));

figure
for ei=3;%1:4%:3;%4;
    exp=experiments{ei};
    
    for binSizei=1:length(binSize_tested);
        binSize=binSize_tested(binSizei);
        
        for lagi=1%:length(lags_tested);
            lags=lags_tested{lagi};
            
            load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/SLeach/regression_SLeach_lag-10-10' ],'sig_fdr');
            sig_betaClass=sig_fdr;
            
            load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
            ris=cell2mat(cellfun(@(x) find(ismember(rnames_table,x)),rnames,'UniformOutput',0));
            r2_sl=nan([length(rnames_table) size(r2,2)]);
            r2_sl(ris,:)=r2;
            
            load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/LLselfother/regression_LL_binSize' num2str(binSize) '_lag0-0' ],'r2','rnames');
            ris=cell2mat(cellfun(@(x) find(ismember(rnames_table,x)),rnames,'UniformOutput',0));
            r2_ll=nan([length(rnames_table) size(r2,2) size(r2,3)]);
            r2_ll(ris,:,:)=r2;
            r2_llm=nanmean(r2_ll,3);
            
            keptT=[min(find(nansum(r2_sl,1)>0)):max(find(nansum(r2_sl)>0))];
            herd(:,1)=corr_col(r2_sl(:,keptT)',r2_llm(:,keptT)');
            
            for perm=1:48;
                if exist([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/perm/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_perm' num2str(perm) '.mat']);
                    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/perm/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_perm' num2str(perm)],'r2','rnames');
                    ris=cell2mat(cellfun(@(x) find(ismember(rnames_table,x)),rnames,'UniformOutput',0));
                    r2_sl=nan([length(rnames_table) size(r2,2)]);
                    r2_sl(ris,:)=r2;
                    
                    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/LLselfother/perm/regression_LL_binSize' num2str(binSize) '_lag0-0' '_perm' num2str(perm)],'r2','rnames');
                    ris=cell2mat(cellfun(@(x) find(ismember(rnames_table,x)),rnames,'UniformOutput',0));
                    r2_ll=nan([length(rnames_table) size(r2,2) size(r2,3)]);
                    r2_ll(ris,:,:)=r2;
                    r2_llm=nanmean(r2_ll,3);
                    
                    herd_null(:,perm)=corr_col(r2_sl(:,keptT)',r2_llm(:,keptT)');
                end
            end
            
            herd_z=real(atanh(herd));
            herd_null_z=real(atanh(herd_null));
            
            for ri=1:length(rnames_table);
                rname=rnames_table(ri);
                %   if    sig_betaClass(ri)==1;
                [~,p(ri,1),~,stats]=ttest(herd_null_z(ri,:),herd_z(ri),'tail','left');
                t(ri,1)=-stats.tstat;
                %   else
                %     p(ri,1)=NaN;
                %    t(ri,1)=NaN;
                %   end
            end
            
            p_adj=nan(size(p));
            [~,~,p_adj(~isnan(p))]=fdr(p(~isnan(p)));
            
            rnames_table(p_adj<.05)
           subplot(3,1, binSizei)
            boxplot(herd_null_z');
            % scatter(repmat(1:61,1,48),herd_null_z(:),20,'k','filled','MarkerFaceAlpha',.1);
            hold on
            scatter(1:61,herd_z,40,'r','filled')
            set(gca,'xtick',1:61,'xticklabels',rnames_table);
            xtickangle(45)
            xlim([0 61]);
            ylabel({exp, ['Bin Size:' num2str(binSize) ]})
            grid on
            text(find(p_adj<.05),herd_z(p_adj<.05),'*','fontsize',14);
            ylim([-0.25 1])
            
            %             save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd/herd_regression_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'herd','herd_null','p_adj','t','p','rnames');
            %
            %             nii=roiTable2wholeBrainNii_mor([roi_ids(p_adj<.05), t(p_adj<.05)]);
            %             save_nii(nii,[expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd/herdt_regression_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_fdr.nii']);
            %
        end
    end
end

