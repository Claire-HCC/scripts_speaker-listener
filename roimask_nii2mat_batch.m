clearvars;
set_parameters;

rfs=cellstr(ls([expdir 'roi_mask\nifti\*R.nii']));
crop_start=0;
voln=[];

for ri = 1:length(rfs);
    rf = [expdir 'roi_mask\nifti\' rfs{ri}];
    
    [data, datasize]=nii2mat(rf,crop_start,voln);
    data=logical(data);
    rf_new=strrep(strrep(rf,'nii','mat'),'nifti','mat');
    save(rf_new,'data','-v7.3');
    clear data datasize;
end


