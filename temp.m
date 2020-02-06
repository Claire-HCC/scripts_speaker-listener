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

for binSizei=1;%1:length(binSize_tested);
    binSize=binSize_tested(binSizei);
    
    for lagi=[1 2 4];%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        figure
        for ei=1;%1:4;
            exp=experiments{ei};
            
            load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
            r2_sl_b=r2;
            ri=find(ismember(rnames,'HG_L'));
            
            
            load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_binSize' num2str(binSize) '_lag0-0' ],'r2');
            r2_ll_b=nanmean(r2,3);
            listenerN=size(r2,3);
            
            keptT=[min(find(nansum(r2_sl_b,1)>0)):max(find(nansum(r2_sl_b)>0))];
            
            subplot(1,2,1)
            scatter(r2_sl_b(ri,keptT)',r2_ll_b(ri,keptT)',40,'filled');
            title('bin')
            
            
            load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2_byTime','rnames');
            r2_sl=r2_byTime;
            ris=cell2mat(cellfun(@(x) find(ismember(rnames_table,x)),rnames,'UniformOutput',0));
            ri=find(ismember(rnames,'dPCUN'));
            
            
            load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag0-0' ],'r2_byTime');
            r2_ll=nanmean(r2_byTime,3);
            listenerN=size(r2,3);
            
            keptT=[min(find(nansum(r2_sl,1)>0)):max(find(nansum(r2_sl)>0))];
            
            subplot(1,2,2)
            scatter(r2_sl(ri,keptT)',r2_ll(ri,keptT)',40,'filled');
            
            
            
        end
        
    end
end

