clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
lags_tested={-4:4,-10:10};


for lagi=1:length(lags_tested);
    lags=lags_tested{lagi};
    
    figure;
    for ei=1:4;
        exp=experiments{ei};
        load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'pfdr');
        load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2');
        load([expdir exp '/bhv/comprehensionScore.mat'],'score');
        % gdata(:,:,subjects_excluded{ei})=NaN;
        r2=squeeze(nanmean(r2(pfdr<.05,:,:),1));
        
        subjs=find(~isnan(score) & ~isnan(r2));
        [r p]=corr(r2(subjs),score(subjs),'type','spearman')
        
        subplot(2,2,ei);
        scatter(r2,score,20,'k','filled');
        title({['lag' num2str(min(lags)) '-' num2str(max(lags))],exp,sprintf('spearman R=%.2f, p=%.3f',r,p)});
        xlabel('Speaker-Listener couplong');
        ylabel('Comprehension score');
        
    end
end



