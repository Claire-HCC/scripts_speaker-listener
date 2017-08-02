loc='mypc';
% loc='cluster';
if strcmp(loc,'mypc');
    expdir='Y:\claire\speaker-listener\';
else
    pcdir='';
end

addpath([expdir '\scripts\' ]);
experiments={'merlin','sherlock'};
tr=[1.5 1.5];
wav_crop_start=[25 25]; %sec