
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

load([expdir '\roi_mask\mor\' 'roi_id_region.mat'],'roi_table');


load([expdir '/bronx/fmri/pattern_regression/' timeUnit '/roi/' froidir '/bz_byHerdCluster.mat'],'b_sl','b_ll','herd_sig','rnames','lags');

i=(sum(herd_sig==0,2)==0);
rnames=rnames(i);
b_sl=b_sl(i,:,:);

% 
% b_sl_z=zscore(b_sl,0,2);
% for ri =1:length(rnames);
%     rname=rnames{ri};
% 
%         
%         figure;
%         plot(lags*1.5,squeeze(b_sl_z(ri,:,1)),'linewidth',2,'color',[0.3 0.7 0.3]);
%         hold on
%         plot(lags*1.5,squeeze(b_sl_z(ri,:,2)),'linewidth',2,'color',[0.7 0.3 0.3]);
%         hold off
%         xlabel(['speaker precedes      shift (sec)    listener precedes']);
%         ylabel('mean weigths (normalized)');
%         set(gca,'fontsize',14)
%         grid on
%         
%         [r(ri) p(ri)]=corr(b_sl_z(ri,:,1)',b_sl_z(ri,:,2)');
%       %  if (p<(0.05/length(rnames))); sig='*'; else sig =''; end
%         
%         title({strrep(rname,'_',' '),['R=' sprintf('%.02f',r(ri)) ]});
%         
%         print(gcf,[expdir '/bronx/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herdShared_' rname '.tif'],'-dtiff');
% 
% end



for ei=1:2;
    exp=experiments{ei};

    b_sl_z=zscore(b_sl(:,:,ei),0,2);
    b_clusters(:,ei)=kmeans_opt(b_sl_z);


%     for ci=1:max(b_clusters);
%         roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(b_clusters==ci));
%         roi_ids=cell2mat(roi_table.id(roi_table_inds));
%         nii=roiTable2wholeBrainNii_mor([roi_ids, ci*ones(length(roi_ids),1)]);
%         save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herdShared_cluster' num2str(ci) '.nii']);
%     end
%
    for ci=1:max(b_clusters);
        bz_byHerdSharedCluster(ci,:,ei)=mean(b_sl_z(b_clusters(:,ei)==ci,:),1);
        figure;
        plot(lags*tr(ei),mean(b_sl_z(b_clusters(:,ei)==ci,:),1),'k','linewidth',2);
        xlabel(['speaker precedes      shift (sec)    listener precedes']);
        ylabel('mean weigths');
        set(gca,'fontsize',14)
        grid on
        title([exp ': cluster' num2str(ci)]);
        print(gcf,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herdShared_Cluster' num2str(ci) '.tif'],'-dtiff');

    end
end
% 
% [r p]=corr(bz_byHerdSharedCluster(:,:,1)',bz_byHerdSharedCluster(:,:,2)')
% figure;
% 
% % bz_byHerdSharedCluster(:,:,2)=bz_byHerdSharedCluster([2 3 1],:,2);
% for ci=1:3;
%     figure;
%     plot(lags*tr(ei),squeeze(bz_byHerdSharedCluster(ci,:,1)),'linewidth',2,'color',[0.3 1 0.3]);
%     hold on
%      plot(lags*tr(ei),squeeze(bz_byHerdSharedCluster(ci,:,2)),'linewidth',2,'color',[1 0.3 0.3]);
%     hold off
%      xlabel(['speaker precedes      shift (sec)    listener precedes']);
%     ylabel('mean weigths');
%     set(gca,'fontsize',14)
%     grid on
% end
% legend(experiments(1:2));
% legend boxoff
% 
% 
% 



