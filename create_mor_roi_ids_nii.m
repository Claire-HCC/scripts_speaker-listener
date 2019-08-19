clear all

set_parameters;
load([expdir '\roi_mask\mor\' 'roi_id_region.mat'],'roi_table');

% roi_ids=zeros(voxn,1);
% 
% for ri=1:length(roi_table.id);
%     rname=roi_table.region{ri};
%     rid=roi_table.id{ri};
%     load([expdir '\roi_mask\mor\mat\' rname '.mat'],'roimask');
%     roi_ids(roimask==1)=rid;
% end
% 
% save([expdir '\roi_mask\mor\roi_ids.mat'],'roi_ids');
% nii=mat2nii(roi_ids);

%save_nii(nii,[expdir '\roi_mask\mor\roi_ids.nii']);


load([expdir '\roi_mask\mor\roi_ids.mat'],'roi_ids');
categories=unique(roi_table.category);
for ci=1:length(unique(roi_table.category));
    category=categories{ci};
    rids_selected=cell2mat(roi_table.id(ismember(roi_table.category,category)));
    
    roimask=zeros(voxn,1);
    roimask(ismember(roi_ids,rids_selected))=1;
    save([expdir '\roi_mask\mor\mat\' category '.mat'],'roimask');
    
    nii=mat2nii(roimask);
    save_nii(nii,[expdir '\roi_mask\mor\nii\' category '.nii']);
end

for ri=2:length(roi_table.id);
    rname=roi_table.region{ri};
    rid=roi_table.id{ri};
    roimask=zeros(voxn,1);
    roimask(roi_ids==rid)=1;
      save([expdir '\roi_mask\mor\mat\' rname '.mat'],'roimask');
    
    nii=mat2nii(roimask);
    save_nii(nii,[expdir '\roi_mask\mor\nii\' rname '.nii']);
end
