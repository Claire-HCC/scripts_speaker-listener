% the "uncropped" speaker data is really uncropped.
% For the "uncropped" listeners's data, actually the 5 scans are already
% cropped. 2 for the delay between first trigger and
% stimuli start. 3 for hrf delay. 

% here I cropped 20 scans from the start for the listeners. That's the
% opening music and silence.
% for the speaker, I found the gap between scan start and the recall based
% on the sound recordings and add 3 scans to crop to account for hrf delay.

% clear all
% close all
% set_parameters
%
% for ei=1:2;
%
%     exp= experiments{ei};
%
%     switch exp
%         case 'merlin';
%             crop_start=20;
%             voln=585;
%         case 'sherlock';
%             crop_start=20;
%             voln=695;
%     end
%
%
%     subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/uncropped/listener*mat']));
%
%     for si=1:length(subjects);
%         subj=subjects{si};
%         load([expdir experiments{ei} '/fmri/mat/wholeBrain/uncropped/' subj]);
%
%         if exist('subj_data');
%             data=subj_data;
%             save([expdir experiments{ei} '/fmri/mat/wholeBrain/uncropped/' subj ],'data');
%         end
%
%         data=data(:,(crop_start+1):(voln+crop_start));
%         save([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj ],'data');
%
%     end
% end


clear all
close all
set_parameters

% for ei=1:2;
% 
%     exp= experiments{ei};
% 
%     switch exp
%         case 'merlin';
%             crop_start=5+3; % 3 for hrf
%             voln=585;
%         case 'sherlock';
%             crop_start=3+3;  % 3 for hrf
%             voln=695;
%     end
% 
% 
%     subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/uncropped/speaker*mat']));
% 
%     for si=1:length(subjects);
%         subj=subjects{si};
%         load([expdir experiments{ei} '/fmri/mat/wholeBrain/uncropped/' subj]);
% 
%         if exist('subj_data');
%             data=subj_data;
%             save([expdir experiments{ei} '/fmri/mat/wholeBrain/uncropped/' subj ],'data');
%         end
% 
%         data=data(:,(crop_start+1):(voln+crop_start));
%         save([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj ],'data');
% 
%     end
% end

