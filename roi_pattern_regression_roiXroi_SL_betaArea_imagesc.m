
% loc='cluster';
clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
[networks,~,networki]=unique(roi_table.network);
netmat=zeros(length(rnames),length(rnames));
for i=1:length(network);
    netmat(networki==i,networki==i)=i;
end
mat_boundary=(bwboundaries(netmat));

crop_start=10;
lags_tested={-10:10, -20:20, -30:30, -10:-4, -20:-4, -30:-4, -10:-1};



figure;
subploti=reshape(1:8,4,2)';
for ei=3;%
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/SLeach/regression_roiXroi_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_bLagArea'  ],'rnames','blag_pos','blag_neg');
        
        blag_neg_mat=nan(length(rnames),length(rnames));
        for i=1:length(network);
            for j=1:length(network);
                blag_neg_mat(networki==i,networki==j)=nanmean(nanmean(blag_neg(networki==i,networki==j)));
            end
        end
        
        figure
        %      subplot(2,4,subploti(1,ei));
        subplot(1,2,1);
        imagesc(blag_neg,[-5 5]);
        hold on
        visboundaries(mat_boundary,'color','r');
        %       subplot(2,4,subploti(2,ei));
        subplot(1,2,2);
        imagesc(blag_neg_mat,[-5 5]);
        colormap jet
        colorbar
        hold on
        visboundaries(mat_boundary,'color','r');
        title(exp)
        ylabel('speaker rois');
        xlabel('listener rois');
    end
end


figure;
for i=1:6;
subplot(2,3,i)
plot(-10:10,squeeze((nanmean(b(11,networki==i,2:end,:),4)))','linewidth',2);
legend(rnames(networki==i));
title(networks{i});
% legend boxoff
ylabel('beta values');
xlabel('speaker precedes---TR---listeners precede')
xlim([-10 20]);
grid on
end



figure;
for i=1:6;
subplot(2,3,i)
plot(-10:10,squeeze((nanmean(r(11,networki==i,:,:),4)))','linewidth',2);
legend(rnames(networki==i));
title(networks{i});
% legend boxoff
ylabel('beta values');
xlabel('speaker precedes---TR---listeners precede')
xlim([-10 20]);
grid on
end









