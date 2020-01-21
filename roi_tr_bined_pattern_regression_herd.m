clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

binSize_tested=[10 15 30]; % tr;
lags_tested={ -10:-4, -20:-4, -30:-4, -10:-1};
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
    rnames_table=table2array(roi_table(:,3));
    
herd=nan([61 4]);
herd_null=nan([61 48 4]);
for binSizei=3;%1:length(binSize_tested);
    binSize=binSize_tested(binSizei);
    
    for lagi=3;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        for ei=1:4;
            exp=experiments{ei};

            load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
            r2_sl=r2;
              ris=cellfun(@(x) find(ismember(rnames_table,x)),rnames);
              
            
            load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_binSize' num2str(binSize) '_lag0-0' ],'r2');
            r2_ll=nanmean(r2,3);
            listenerN=size(r2,3);
            
            keptT=[min(find(nansum(r2_sl,1)>0)):max(find(nansum(r2_sl)>0))];
            herd(ris,ei)=corr_col(r2_sl(:,keptT)',r2_ll(:,keptT)');
            
            for perm=1:listenerN;
                load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/perm/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  '_perm' num2str(perm)],'r2');
                r2_s_temp=r2;
                
                load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/perm/regression_LL_binSize' num2str(binSize) '_lag0-0' '_perm' num2str(perm)],'r2');
                r2_l_temp=r2;
                
               herd_null(ris,perm,ei)=corr_col(r2_s_temp(:,keptT)',r2_l_temp(:,keptT)');
            end
            
            herd_z=real(atanh(herd(:,ei)));
            herd_null_z=real(atanh(herd_null(:,:,ei)));
            
            for ri=1:length(rnames);
                [~,p(ri,1),~,stats]=ttest(herd_null_z(ri,:),herd_z(ri));
                t(ri,1)=-stats.tstat;
            end
            % p=sum(herd_null>herd,2)/listenerN;
            
            subplot(4,1,ei)
            % boxplot(herd_null_z');
            scatter(repmat(1:61,1,48),herd_null_z(:),20,'k','filled','MarkerFaceAlpha',.2);
            hold on
            scatter(1:61,herd_z,40,'r','filled')
            set(gca,'xtick',1:61,'xticklabels',rnames_table);
            xtickangle(45)
            xlim([0 61]);
            ylabel(exp)

        end
        
    end
end

