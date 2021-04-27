clear all
close all
set_parameters;
for ei=4;
    exp=experiments{ei};
    load([expdir exp '/fmri/timeseries/tr/wholeBrain/3x3x3mm/speaker.mat'],'data');
    nii=mat2nii_3x3x3mm(data);
    save_nii(nii,[expdir exp '/fmri/timeseries/tr/wholeBrain/3x3x3mm/speaker01.nii']);
    
    load([expdir exp '/fmri/timeseries/tr/wholeBrain/3x3x3mm/listenerAll.mat'],'data');
    for si=1:size(data,3);
        nii=mat2nii_3x3x3mm(data(:,:,si));
        save_nii(nii,sprintf('%s%s/fmri/timeseries/tr/wholeBrain/3x3x3mm/listener%02d.nii',expdir,exp,si));
    end
end

