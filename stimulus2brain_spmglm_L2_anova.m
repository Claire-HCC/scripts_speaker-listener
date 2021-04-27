% spm12
clear all
% loc='cluster';
set_parameters
modelNames={'model_log','model_linear_noVolterra'};
conds={'syllable','word','sentence','paragraph','wd_sylN','sn_wdN','pr_snN'};

% addpath(spm_path); % comment to run in compiled spm8
spm('defaults', 'FMRI');
spm_jobman('initcfg');

for mi=2;%1:length(modelNames);
    
    for condi=4;%:length(conds);
        modelName=modelNames{mi};
        cond=conds{condi};
        modeldir=sprintf('%s/sherlock/fmri/temporal/stimulus2brain_glm/%s/SherlockMerlin/%s',expdir,modelName,cond);
        mkdir(modeldir);
        matlabbatch{1}.spm.stats.factorial_design.dir = {modeldir};
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'story';
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Subjects';
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
        %%
        load(sprintf('%s/sherlock/fmri/temporal/stimulus2brain_glm/%s/L1_spm/s01/SPM.mat',expdir,modelName));
        cnames={SPM.xCon.name};
        contri=find(ismember(cnames,cond));
        
        clear f
        fi=1;
        for ei=1:2;
            exp=exp_parameters.experiments{ei};
            for si=1:18;
                f{fi,1}=sprintf('%s%s/fmri/temporal/stimulus2brain_glm/%s/L1_spm/s%02d/scon_%04d.nii,1',expdir,exp,modelName,si,contri);
                fi=fi+1;
            end
        end
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans =f;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.imatrix = [1 1 1 1
            2 1 2 1
            3 1 3 1
            4 1 4 1
            5 1 5 1
            6 1 6 1
            7 1 7 1
            8 1 8 1
            9 1 9 1
            10 1 10 1
            11 1 11 1
            12 1 12 1
            13 1 13 1
            14 1 14 1
            15 1 15 1
            16 1 16 1
            17 1 17 1
            18 1 18 1
            19 2 1 1
            20 2 2 1
            21 2 3 1
            22 2 4 1
            23 2 5 1
            24 2 6 1
            25 2 7 1
            26 2 8 1
            27 2 9 1
            28 2 10 1
            29 2 11 1
            30 2 12 1
            31 2 13 1
            32 2 14 1
            33 2 15 1
            34 2 16 1
            35 2 17 1
            36 2 18 1];
        %%
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum = 2;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em ={sprintf('%s/sherlock/fmri/temporal/stimulus2brain_glm/%s/SherlockMerlin/L1_masks.nii',expdir,modelName)};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'pos';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 1 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111 0.111111111111111];
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'pos_sherlock';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 0 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556];
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'pos_merlin';
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 1 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556 0.0555555555555556];
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.delete = 0;
         save([modeldir '/matlabbatch.mat'],'matlabbatch')
      spm_jobman('run',matlabbatch);
    end
end