

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment-specific parameters
clearvars '*'
loc='mypc';
set_parameters
modelName= 'event_anova';
roles={'listener','speaker01'}; % 'speaker'
exp='merlin';
modeldir=[expdir exp '/fmri/nii/wholeBrain/group/' modelName ];

% addpath(spm_dir); % comment to run in compiled spm8
spm('defaults', 'FMRI');
spm_jobman('initcfg');
clear matlabbatch
%%%%%%%%%%%%%%%%

matlabbatch{1}.spm.stats.con.spmmat = {[modeldir '/SPM.mat']};
matlabbatch{1}.spm.stats.con.delete = 1;

n=1;
for evi=1:22;
    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = ['event' num2str(evi)]
    matlabbatch{1}.spm.stats.con.consess{n}.tcon.convec =  [-1/21*ones(1,evi-1) 1 -1/21*ones(1,22-evi) ]
    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
    n=n+1;
end

matlabbatch{1}.spm.stats.con.consess{end+1}.fcon.name = 'allF';
matlabbatch{1}.spm.stats.con.consess{end}.fcon.convec = {[eye(22) 1/18*ones(22,18) ]}';
matlabbatch{1}.spm.stats.con.consess{end}.fcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{end+1}.fcon.name = 'allDiff';
matlabbatch{1}.spm.stats.con.consess{end}.fcon.convec = {[eye(21) zeros(21,1)]-[zeros(21,1) eye(21) ]  }';
matlabbatch{1}.spm.stats.con.consess{end}.fcon.sessrep = 'none';

% matfile = sprintf('%sContra_%s_%s_%s.mat', batch_dir,subjects{csubj},condition,modelName);
% save(matfile,'matlabbatch');
% jobs{csubj} = matfile;
spm_jobman('run',matlabbatch);
clear matlabbatch








