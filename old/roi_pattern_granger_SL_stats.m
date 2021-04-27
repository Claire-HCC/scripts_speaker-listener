clear all
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
crop_start=10;
lags_tested={-10:-4, -20:-4, -30:-4, -10:-1};

figure
for ei=[1:4]%1:2;
    exp=experiments{ei};
    listenerN=listenerNs(ei);
    
    rnames=table2array(roi_table(:,3));
    ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
    rnames=rnames(ris);
    
    for lagi=3;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/granger_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'F_s2l','rnames');
        F_s2l_real=F_s2l;
        
        for perm=1:listenerN;
            load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/perm/granger_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm' num2str(perm)],'F_s2l','rnames');
            F_s2l_null(:,perm)=F_s2l;
        end
        %  clear F_s2l_null
        p=sum(F_s2l_null>F_s2l_real,2)/listenerN;
        rnames(p<0.05);
subplot(4,1,ei)
        boxplot(F_s2l_null'); hold on
        scatter(1:length(rnames),F_s2l_real,40,'k','filled'); hold off
        title([exp '; lag:' num2str(min(lags)) '-' num2str(max(lags))])
    end
    clear F_s2l_null
end
