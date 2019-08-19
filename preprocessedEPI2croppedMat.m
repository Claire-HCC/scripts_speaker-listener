clearvars;
set_parameters;

roles={'speaker','listener'};
timeTable=readtable([expdir 'scripts_speaker-listener/crop_start_voln.csv']);

tic % 26 subjects take approximately 30 min with tr as timeUnit.
% for ri=1%:length(roles);
%     role=roles{ri};
%     for ei=3:4;%3;
%         exp=experiments{ei};
%         subjects=cellstr(ls([expdir '/fmri_preprocesing/subjects/' exp '_' role '*' ]));
%         
%         
%         for si = 1%1:length(subjects);
%             subj=subjects{si};
%             % ther eare 3 versions of pieman2, we chose version 3 in the end
%             if strcmp(exp,'pieman2') & strcmp(role,'speaker'); ses=3; else ses=1;end
%             
%             fnii = sprintf('%s\\fmri_preprocesing\\subjects\\%s\\analysis\\preproc\\preproc%02d.feat\\trans_filtered_func_data.nii',expdir,subj,ses);
%             
%             ti=find(strcmp(experiments{ei},timeTable.experiment) & strcmp(strrep(subj,[exp '_'],''),timeTable.subject));
%             % no crop
%             %  [data, datasize]=nii2mat(fnii);
%             
%             [data, datasize]=nii2mat(fnii,timeTable.crop_start(ti),timeTable.voln(ti));
%             fmat = sprintf('%s\\%s\\fmri\\mat\\wholeBrain\\%s%02d.mat',expdir,experiments{ei},role,si);
%             save(fmat,'data','datasize','-v7.3');
%             
%         end
%         
%     end
% end

for ri=1%:2;%1:length(roles);
    role=roles{ri};
    
    for ei=1%:2;
        exp= experiments{ei};
        subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/uncropped/' role '*mat']));
        
        for si=1:length(subjects);
            subj=strrep(subjects{si},'.mat','');
            load([expdir experiments{ei} '/fmri/mat/wholeBrain/uncropped/' subj '.mat']);
            
            ti=find(strcmp(experiments{ei},timeTable.experiment) & strcmp(subj,timeTable.subject));
            s=(timeTable.crop_start(ti)+1);
            e=s+timeTable.voln(ti)-1;
            
            
            %         if exist('subj_data');
            %             data=subj_data;
            %             save([expdir experiments{ei} '/fmri/mat/wholeBrain/uncropped/' subj ],'data');
            %         end
            
            data=data(:,s:e);
            save([expdir experiments{ei} '/fmri/mat/wholeBrain/' subj '.mat'],'data');
            
        end
    end
end

