clear all
close all
k=6;

set_parameters;

for ei=1%:2;
    exp=experiments{ei};
    
    load([expdir exp '/fmri/mat/wholeBrain/isfc/isfc_slm_mean.mat' ],'isfc','vox_selected');
    vox_selected=find(vox_selected);
    figure
    subplot(1,4,1);
    imagesc(isfc,[-0.23 0.23]);
    colormap('jet');
    colorbar;
    
    set(gca,'XAxisLocation', 'top');
    xlabel('speaker');
    ylabel('listener');
    set(gca, 'XTick', []); % center x-axis ticks on bins
    set(gca, 'YTick', []); % center y-axis ticks on bins
    
    % kmeans cluster row vectors. isfc is a listener x speaker corr matrix
    [idx_s,C_s,sumd_s, D_s]  = kmeans(isfc',k);
    [idx_l,C_l,sumd_l, D_l]  = kmeans(isfc,k);
    
    [idx i_s]=sort(idx_s);
    [idx i_l]=sort(idx_l);
    isfc_clustered1=isfc(i_l,:);
    isfc_clustered2=isfc(:,i_s);
    isfc_clustered3=isfc(i_l,i_s);
    
    subplot(1,4,2);
    imagesc(isfc_clustered1,[-0.23 0.23]);
    colormap('jet');
    colorbar;
    set(gca,'XAxisLocation', 'top');
    xlabel('speaker');
    ylabel('listener');
    set(gca, 'XTick', []); % center x-axis ticks on bins
    set(gca, 'YTick', []); % center y-axis ticks on bins
    
    subplot(1,4,3);
    imagesc(isfc_clustered2,[-0.23 0.23]);
    colormap('jet');
    colorbar;
    set(gca,'XAxisLocation', 'top');
    xlabel('speaker');
    ylabel('listener');
    set(gca, 'XTick', []); % center x-axis ticks on bins
    set(gca, 'YTick', []); % center y-axis ticks on bins
    % xticklabels(1:k);
    % yticklabels(1:k);
    
    subplot(1,4,4);
    imagesc(isfc_clustered3,[-0.23 0.23]);
    colormap('jet');
    colorbar;
    set(gca,'XAxisLocation', 'top');
    xlabel('speaker');
    ylabel('listener');
    set(gca, 'XTick', []); % center x-axis ticks on bins
    set(gca, 'YTick', []); % center y-axis ticks on bins
    
    % better show with xjview
    for ki=1:k;
        cluster_s=zeros(voxn,1);
        cluster_l=zeros(voxn,1);
        % vi_1=vi(sort(idx2));
        
        cluster_s(vox_selected(idx_s==ki))=ki;
        cluster_l(vox_selected(idx_l==ki))=ki;
        nii=mat2nii(cluster_s);
        save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/isfc/' exp '_isfc_cluster_s_cluster' num2str(ki) '.nii']);
        
        nii=mat2nii(cluster_l);
        save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/isfc/' exp  '_isfc_cluster_l_cluster' num2str(ki) '.nii']);
    end
    
end