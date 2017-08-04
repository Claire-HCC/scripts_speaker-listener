loc='mypc';
%loc='cluster';
if strcmp(loc,'mypc');
    expdir='Y:\claire\speaker-listener\';
else
    % expdir='/scratch/claire/speaker-listener/';
    expdir='/scratch/claire/speaker-listener/';
end

addpath([expdir '\scripts_speaker-listener\' ]);
experiments={'merlin','sherlock'};
tr=[1.5 1.5];
wav_crop_start=[25 25]; %sec