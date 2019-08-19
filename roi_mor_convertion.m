%  for f in $fs 
% do
% f_new=`echo $(basename $f)`
% flirt -in $f -ref ../MNI152_T1_3mm_brain.nii -out ./nii/$f_new -applyxfm
% fslmaths ./nii/$f_new -thr 0.9 -bin ./nii/$f_new
% fslchfiletype NIFTI ./nii/$f_new
% done
table=readtable('../roi_id_region.txt');
fs=cellstr(ls('*nii'));

for fi=1:61;
f=fs{fi};
id=regexp(f,'(roi[0-9]*)_','tokens');
ri=find(strcmp(table.id,id{1}{1}));
roi_name=table.region(ri);
copyfile(f,[roi_name{1} '.nii']);

nii=load_nii(f);
roimask=nii.img(:);
save(['../mat/' roi_name{1} '.mat'],'roimask');
end