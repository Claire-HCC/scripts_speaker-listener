clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10, -10:-1, 0, 1:10};

figure;
for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) '_comprehensionscore' ],'r2','score','r','p','pfdr');
        subplot(2,2,ei);
        histogram(r,[-1:0.1:1]);
        hold on;
rnames(pfdr<.05)
    end
    title(exp)
end
legend(cellfun(@(x) sprintf('lags %d-%d',min(x),max(x)),lags_tested,'Uniformoutput',0))
legend boxoff



%
%   scatter(r2,score);
%
%         c = polyfit(r2,score,1);
%         y = polyval(c,r2);
%         plot(r2,y);
%         title(sprintf('',exp'r,p))