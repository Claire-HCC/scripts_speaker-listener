clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10, -10:-1, 0, 1:10};

% for ei=1:4;
%     exp=experiments{ei};
%     
%     for lagi=1:length(lags_tested);
%         lags=lags_tested{lagi};
%         
%         load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags))],'r2');
%         load([expdir exp '/bhv/comprehensionScore.mat'],'score');
%         r=nan(size(rnames));
%         p=nan(size(rnames));
%         
%         for ri=1:length(rnames);
%             subjs=find(~isnan(score) & ~isnan(squeeze(r2(ri,:))'));
%             if ~isempty(subjs);
%                 [r(ri,1), p(ri,1)]=corr(score(subjs),squeeze(r2(ri,subjs)'),'type','spearman','tail','right');
%             end
%         end
%         pfdr=nan(size(p));
%         [~,~,pfdr(~isnan(p))]=fdr(p(~isnan(p)));
%         
%         save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) '_comprehensionscore' ],'r2','score','r','p','pfdr');
%         
%         nii=roiTable2wholeBrainNii_mor([roi_ids, r]);
%         save_nii(nii,[expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/'  froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags))  '_comprehensionScore.nii']);
%         
%     end
% end


ri=18;
for lagi=3;%1:length(lags_tested);
    lags=lags_tested{lagi};
    
    figure;
    for ei=1:4;
        exp=experiments{ei};
        load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) '_comprehensionscore' ],'r2','score','r','p','pfdr');
        
        subplot(2,2,ei);
        scatter(r2(ri,:),score,40,'filled','k');
        title({[exp ' ,' rnames{ri}],sprintf('R=%.2f, p=%.3f',r(ri),p(ri))});
        xlabel({'SL coupling', '(r-squared of pattern regression model', sprintf('with lag%d-%d)',min(lags),max(lags))});
        ylabel('Comprehension Score');
        set(gca,'fontsize',12);
        
    end
end


