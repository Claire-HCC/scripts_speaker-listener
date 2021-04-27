% spm12
clear all
set_parameters
modelName='model_log_withoutPause';
% addpath(spm_path); % comment to run in compiled spm8
spm('defaults', 'FMRI');
spm_jobman('initcfg');

for ei=1;
    exp=experiments{ei};
    load(sprintf('%s/%s/fmri/temporal/stimulus2brain_glm/%s/L1/L1.mat',expdir,exp,modelName),'xnames');
    
    for bi=1:length(xnames);
        xname=xnames{bi};
        mkdir([expdir '\'  exp '\fmri\temporal\stimulus2brain_glm\' modelName '\' xname]);
        matlabbatch{1}.spm.stats.factorial_design.dir = {[expdir '\'  exp '\fmri\temporal\stimulus2brain_glm\' modelName '\' xname]};
        %%
        for si=1:18;
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{si,1} = sprintf('%s/%s/fmri/temporal/stimulus2brain_glm/%s/L1//%s_subj%02d.nii',expdir,exp,modelName,xname,si);
        end
        %%
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
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