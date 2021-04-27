
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
lags_tested={-10:10, -20:20, -30:30, -10:-4, -20:-4, -30:-4, -10:-1};


figure
for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/SLeach/regression_SLeach_lag' num2str(min(lags)) '-' num2str(max(lags)) '_classification'  ],'rnames','acc','sig_fdr','b_s','b_l','p','p_adj','foldN');
        acc_m=nanmean(acc,2);
        
        subplot(4,1,ei)
        boxplot(acc',1:length(rnames));
        set(gca,'xtick',1:61,'xticklabels',rnames);
        xtickangle(45)
        xlim([0 61]);
        ylim([0 1]);
        line([0 61],[0.5 0.5],'color','k')
        title(exp);
        ylabel('Beta profile lassification Acc.')
        grid on
        % text(find(p_adj<.05),max(herd_null_z(p_adj<.05,:),[],2)+0.1,'*','fontsize',14);
        text(find(p_adj<.05),acc_m(p_adj<.05),'*','fontsize',14);
    end
end

