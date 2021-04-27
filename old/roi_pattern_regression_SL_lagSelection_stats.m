clear all
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
crop_start=10;
perms=100;

for ei=[1 2 3 4];%
    exp=experiments{ei};
    rnames=table2array(roi_table(:,3));
    ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
    rnames=rnames(ris);
    
    for perm=1:perms;
        
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/perm/regression_SL_lagSelection_perm' num2str(perm)  ],'lags','rnames','keptT','r2_test','r2_train');
        r2_test_all(:,:,perm)=r2_test;
        r2_train_all(:,:,perm)=r2_train;
        
        % clear r2_train r2_test
    end
    
    %     p=sum(r2_test_all>0,3)/perms;
    r2_test_m=mean(r2_test_all,3);
    
    for ri=1:size(r2_test_all,1);
        [~,p(ri,:),~,stats]=ttest(squeeze(r2_test_all(ri,:,:))',0,'tail','right');
        t(ri,:)=stats.tstat;
        
        if sum(p(ri,:)<(0.05)) >0;
            [tmx,~]=max(t(ri,p(ri,:)<0.05));
            lags_peak(ri,1)=-find(t(ri,:)==tmx);
        else;
            lags_peak(ri,1)=NaN;
        end
    end
    
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lagSelection_stats'  ],'lags','rnames','r2_test_all','lags_peak','t');
    
    subplot(1,4,ei)
    hist(lags_peak);
    title(exp)
    
    clear r2_test_all r2_train_all
end
%  hist(lags_peak(~isnan(lags_peak)),length(lags))
% plot(mean(r2_test_all(ri,:,:),3))
% plot(squeeze(r2_test_all(ri,:,:)))
% imagesc(zscore(t,0,2))
% plot(zscore(t,0,2)'); grid on