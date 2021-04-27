%
% % %% boxplot
% clear all;
%
% % loc='cluster';
% set_parameters;
% timeUnit='tr' ;
% froidir='mor';
%
% binSize_tested=[10 20 30 40]; % tr;
% lags_tested={-10:10 -10:-4, -10:-1, -20:-4, -30:-4,-10:-3};
% load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
% rnames_table=table2array(roi_table(:,3));
% roi_ids=cell2mat(table2array(roi_table(:,1)));
%
% figure
% for ei=3;
%     exp=experiments{ei};
%
%     for lagi=2;%1:length(lags_tested);
%         lags=lags_tested{lagi};
%
%         for binSizei=1:length(binSize_tested);
%             binSize=binSize_tested(binSizei);
%
%             load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/herd/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'herd','herd_null','pfdr','t','p','rnames');
%             herd_z=atanh(herd);
%             herd_null_z=atanh(herd_null);
%
%             rnames_table(pfdr<.05)
%             subplot(2,2, binSizei)
%             boxplot(herd_null_z');
%
%             hold on
%             scatter(1:61,herd_z,40,'r','filled')
%             set(gca,'xtick',1:61,'xticklabels',rnames_table);
%             xtickangle(45)
%             xlim([0 61]);
%             ylabel({exp, ['Bin Size:' num2str(binSize) ],['lag' num2str(min(lags)) '-' num2str(max(lags)) ]})
%             grid on
%             text(find(pfdr<.05),herd_z(pfdr<.05),'*','fontsize',14);
%             ylim([-0.25 1.25])
%
%
%         end
%     end
% end
%
%
%
%
