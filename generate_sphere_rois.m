clear all

loc='mypc';
set_parameters;
rois={'HG_L','Angular_L','precuneus'};

for ei=1:11;
    exp=exp_parameters.experiments{ei};
    
    isc=load_nii([expdir exp '\fmri\temporal\circularlagcorr\tr\vox\LL_leave1out\isc_r.nii']);
    
    % maximum is specified as the centre of the sphere in mm in MNI space
    for ri=1;%1:length(rois);
        roi=rois{ri};
        
        roimask=load_nii([expdir '\roi_mask\aal\nii\' roi '.nii']);
      
        roi_isc=isc.img;
        roi_isc(roimask.img==0)=NaN;
        cor=[];
        [cor(1) cor(2) cor(3)]=ind2sub(volsize,find(roi_isc(:)==max(roi_isc(:))));
        sphere_centre=cor2mni(cor);
        % sphere_centre = cell2mat(table2array(exp_parameters(ei,['isc_peak_' roi])));
        sphere_radius =5; % mm
        sphere_roi = maroi_sphere(struct('centre', sphere_centre,...
            'radius', sphere_radius));
        
        % Save as image
        space=spm_vol([expdir '/roi_mask/MNI152NLin2009cAsym_3x3x4mm_brain.nii']);
        save_as_image(sphere_roi, [expdir '/roi_mask/isc_peak/nii/' roi '_' exp '.nii'],space);
        
        nii=load_nii( [expdir '/roi_mask/isc_peak/nii/' roi '_' exp '.nii']);
        roimask=nii.img(:);
        save([expdir '/roi_mask/isc_peak/mat/' roi '_' exp '.mat'],'roimask');
    end
end
