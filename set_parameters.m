
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
% addpath([expdir '/scripts_speaker-listener/FastICA_25/' ]);


volsize=[65 77 49];
voxn=prod(volsize);

load([expdir '/scripts_speaker-listener/exp_parameters.mat']);
% exluded for comprehensino score < 1.5 SD below the mean (pieman, bronx, black, forgot)
% piemen subj 22 is excluded for translation along z-axis > 3mm
% Asieh has removeed merlin/sherlock subjects with movements > 3mm
% I don't have the movement parameteters of prettymouth and pieman0

% these numbers were chosen based on the audio-envelope analysos
% merlin speaker: 5 tr speaker start delay. 3 Tr hrf correction
% sherlock speaker: 3 tr stimulus start delay and 3 Tr hrf correction
% merlin & sherlock listener: 20 TR (20 music) of the opeining music, (Asieh already cropped 5 TR (2 For stimulus delay (gap between scan start and stimulus start); 3 For hrf correction)); 
% Bronx: crop 21 TRs (8 TR epi/sound gap + 10 silence + 3 hrf). Voiced period lasts for 358 TRs
% Pieman: cropt 17 (8 TR epi/sound gap + 6 TRs of silence + 3 hrf). Voiced period lasts for 267 TRs.
% pieman old: crop 14 TR (10 music + 3 hrf + 1 ??)), voiced period 276 TR
% pieman oldWord (I do not have the sound recording of this, so I just copy these numbes from piemen old, the intact version): crop 14 TR (10 music + 3 hrf + 1 ??)), voiced period 276 TR.
% prettymouth: crop 19 TR (14 music + 3 hrf + 2 ??)). Voiced period lasts for 450 TRs.

