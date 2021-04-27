clear all
% loc='cluster';
set_parameters;
load([expdir '\roi_mask\shen\' 'roi_id_region.mat'],'roi_table');




load([expdir '/roi_mask/shen/roi_ids.mat'],'roi_ids');
networks=unique(roi_table.network);
for ci=1:length(unique(roi_table.network));
    network=networks{ci};
    rids_selected=cell2mat(roi_table.id(ismember(roi_table.network,network)));
    
    roimask=zeros(voxn,1);
    roimask(ismember(roi_ids,rids_selected))=1;
    save([expdir '/roi_mask/shen/mat/network' network '.mat'],'roimask');
    
    nii=mat2nii(roimask);
    save_nii(nii,sprintf('%s/roi_mask/shen/nii/network%s.nii',expdir,network));
end


for ri=1:length(roi_table.id);
    rname=roi_table.region{ri};
    rid=roi_table.id{ri};
    roimask=zeros(voxn,1);
    roimask(roi_ids==rid)=1;
      save([expdir '/roi_mask/shen/mat/' rname '.mat'],'roimask');
    
    nii=mat2nii(roimask);
    save_nii(nii,[expdir '/roi_mask/shen/nii/' rname '.nii']);
end
