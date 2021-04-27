function stimulus2brain_glm_L1(si)
%-----------------------------------------------------------------------
% Job saved on 28-Nov-2020 03:59:35 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6470)
% cfg_basicio BasicIO - Unknown
% -----------------------------------------------------------------------
loc='cluster';
addpath('/jukebox/hasson/claire/spm/spm12');

spm('defaults', 'FMRI');
spm_jobman('initcfg');

set_parameters;
modelName='model_linear_noVolterra';
for ei=1:2;
    exp=exp_parameters.experiments{ei};
    
    modeldir=sprintf('%s/%s/fmri/temporal/stimulus2brain_glm/%s/L1_spm/s%02d',expdir,exp,modelName,si);
    mkdir(modeldir);
    matlabbatch{1}.spm.stats.fmri_spec.dir = {modeldir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.5;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 27;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 14;
    %%
    f=cell(exp_parameters.voln_listener{ei},1);
    for voli=1:exp_parameters.voln_listener{ei};
        f{voli}=sprintf('%s%s/fmri/timeseries/tr/wholeBrain/listener%02d.nii,%d',expdir,exp,si,voli);
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = f;
    
    %%
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {sprintf('%s/%s/fmri/temporal/stimulus2brain_glm/%s/L1_spm/conditionSpec.mat',expdir,exp,modelName)};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};%{sprintf('%s/roi_mask/MNI152NLin2009cAsym_3x3x4mm_brain_gm.nii',expdir)};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'syllable';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'wd_sylN';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'word';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'sn_wdN';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 0 0  1];
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'sentence';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 0  1];
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'pr_snN';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 0 0 1];
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'paragraph';
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 0  0 0 1];
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'noPrLabel';
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [0 0 0 0 0  0 0 1];
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 0;
   save([modeldir '/matlabbatch.mat'],'matlabbatch')
    spm_jobman('run',matlabbatch);
end
