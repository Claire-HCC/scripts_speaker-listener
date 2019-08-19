function nii=roiTable2wholeBrainNii_mor(roiTable,l);
% roiTable: first column is roi labels. Second is the velue in that roi;
% loc='cluster';
set_parameters;
nii=load_nii([expdir '/roi_mask/mor/roi_ids.nii']);
roi_labels_img=nii.img(:);

img=nan(size(roi_labels_img));
for ri=1:length(roiTable(:,1));
    if ~(roiTable(ri,1)==0);
    img(roi_labels_img==roiTable(ri,1))=roiTable(ri,2);
    end
end

nii=mat2nii(img);

end
