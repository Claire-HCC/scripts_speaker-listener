clear all
close all
k=6;

set_parameters;
for ei=1%:2;
    exp=experiments{ei};
    
    load([expdir exp '/fmri/mat/roi/' exp '_isfc_roi.mat' ],'isfc');
    
    figure;
    
    lags=(1:size(isfc,4))-ceil(size(isfc,4)/2);
    figure
    for lagi=1:size(isfc,4);
        isfc_temp=nanmean(isfc(:,:,:,lagi),3);
        
        subplot(ceil(sqrt(size(isfc,4))),ceil(sqrt(size(isfc,4))),lagi);
        imagesc(isfc_temp,[-0.2 0.2]);
        colormap('jet');
        title(['lag ' num2str(lags(lagi)) ' scan']);
        
        set(gca, 'XTick', []); % center x-axis ticks on bins
        set(gca, 'YTick', []); % center y-axis ticks on bins
        
        %     % kmeans cluster row vectors. isfc is a listener x speaker corr matrix
        %     [idx_s,C_s,sumd_s, D_s]  = kmeans(isfc_temp',k);
        %     [idx_l,C_l,sumd_l, D_l]  = kmeans(isfc_temp,k);
        %
        %     [idx i_s]=sort(idx_s);
        %     [idx i_l]=sort(idx_l);
        %     isfc_temp_clustered=isfc_temp(i_l,i_s);
        %
        %     figure
        %     imagesc(isfc_temp_clustered,[-0.2 0.2]);
        %     colormap('jet');
        %     colorbar;
        %
        %     set(gca,'XAxisLocation', 'top');
        %     xlabel('speaker');
        %     ylabel('listener');
        %     set(gca, 'XTick', []); % center x-axis ticks on bins
        %     set(gca, 'YTick', []); % center y-axis ticks on bins
        % xticklabels(1:k);
        % yticklabels(1:k);
        
        %     cluster_s=zeros(size(isfc,1),1);
        %     cluster_l=zeros(size(isfc,1),1);
        %     for ki=1:k;
        %         % vi_1=vi(sort(idx2));
        %
        %         cluster_s(idx_s==ki)=ki;
        %         cluster_l(idx_l==ki)=ki;
        %         %         nii=mat2nii(cluster_s);
        %         %         save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/' exp '_isfc_cluster_s_cluster' num2str(ki) '.nii']);
        %         %
        %         %         nii=mat2nii(cluster_l);
        %         %         save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/' exp  '_isfc_cluster_l_cluster' num2str(ki) '.nii']);
        %     end
    end
    subplot(ceil(sqrt(size(isfc,4))),ceil(sqrt(size(isfc,4))),lagi+1);
    colorbar;
    set(gca,'XAxisLocation', 'top');
    xlabel('speaker');
    ylabel('listener');
end