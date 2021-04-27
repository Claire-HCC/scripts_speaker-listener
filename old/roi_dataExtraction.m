clear all
close all
loc='mypc';
set_parameters

rnames=strrep(cellstr(ls([expdir '/roi_mask/mat/*.mat'])),'.mat','');
roles={'listenerMean','speaker01','listener01','listener01_others'};

for ei=1%:4;
    exp=experiments{ei};
    clear data
    
    for rolei=1:length(roles);
        role=roles{rolei};
        %  subjects=cellstr(ls([expdir experiments{ei} '/fmri/mat/wholeBrain/*'  role '*mat']));
        %    subjects=dir([expdir experiments{ei} '/fmri/mat/wholeBrain/*'  role '*mat']);
        %   subjects={subjects.name};
        %    subjects= subjects(~cellfun(@isempty,(regexp(subjects,'[0-9]*'))));
        
        %     for si =1:length(subjects);
        %    subj=subjects{si};
        epi=load([expdir experiments{ei} '/fmri/mat/wholeBrain/' role '.mat']);
        epi=epi.data;
        epi=zscore(epi,0,2);
        
        data=table;
        for ri=1:length(rnames);
            rname=rnames{ri};
            rmask=load([expdir '/roi_mask/mat/' rname '.mat']);
            rmask=rmask.roimask;
            data.(rname)=nanmean(epi(rmask==1,:))';
        end
        save([expdir experiments{ei} '/fmri/mat/roi/' exp '_' role '_rois.mat' ],'data');
    end
    
end
% end