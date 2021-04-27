% spm12
clear all
set_parameters
modelName='model_log';
% addpath(spm_path); % comment to run in compiled spm8
spm('defaults', 'FMRI');
spm_jobman('initcfg');

for ei=2;
    exp=experiments{ei};
    modeldir=sprintf('%s/%s/fmri/temporal/stimulus2brain_glm/%s/',expdir,exp,modelName);
    load([modeldir  '/L1_spm/s01/SPM.mat']);
    cnames= {SPM.xCon.name};
    
    for ci=1:length(cnames);
                
        clear matlabbatch
        
        cname=cnames{ci};
        mkdir([modeldir '/' cname]);
        
        f=cell(18,1);
        for si=1:18;
            f{si,1}=sprintf('%s/L1_spm/s%02d/con_%04d.nii',modeldir,si,ci);
        end
        
        matlabbatch{1}.spm.stats.factorial_design.dir = {[modeldir '\' cname]};
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = f;
        
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {sprintf('%s/roi_mask/MNI152NLin2009cAsym_3x3x4mm_brain_gm.nii',expdir)};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('modelName estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'pos';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'neg';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.delete = 0;
        
        spm_jobman('run',matlabbatch);

    end
end