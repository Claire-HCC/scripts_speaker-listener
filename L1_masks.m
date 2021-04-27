% spm12
clear all
% loc='cluster';
set_parameters
modelNames={'model_log','model_linear_noVolterra'};

for mi=2;%1:length(modelNames);
    
    modelName=modelNames{mi};
    masks=zeros(voxn,1);
    for ei=1:2;
        exp=exp_parameters.experiments{ei};
        for si=1:18;
            f=sprintf('%s%s/fmri/temporal/stimulus2brain_glm/%s/L1_spm/s%02d/mask.nii',expdir,exp,modelName,si);
            mask=load_nii(f);
            masks=masks+double(mask.img(:));
        end
    end
    
    masks=(masks==(2*18));
    nii=mat2nii(masks);
    save_nii(nii,sprintf('%s/sherlock/fmri/temporal/stimulus2brain_glm/%s/SherlockMerlin/L1_masks.nii',expdir,modelName));
end