
if ~exist('loc');
    loc='mypc';
    % loc='cluster';
end

if strcmp(loc,'mypc');
    expdir='Y:/claire/speaker-listener/';
else
   
    expdir='/scratch/claire/speaker-listener/';
    addpath(genpath([expdir '/../tools/isc_nifti_kit_080916']));
end

addpath([expdir '/scripts_speaker-listener/' ]);
experiments={'pieman','bronx','merlin','sherlock'};
tr=[1.5 1.5 1.5 1.5];

volsize=[65 77 49];
voxn=prod(volsize);

listenerNs=[47 48 18 18]; % listener 10 in pieman did not go through the whole scanning (only 287TRs were obtained)
listeners={[1:9 11:48],[1:48],1:18,1:18};
