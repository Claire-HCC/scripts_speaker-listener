
clear all
close all
set_parameters

for ei=1%:2;
    subjects=dir([expdir experiments{ei} '/fmri/mat/wholeBrain/listener*mat']);
    subjects={subjects.name}';
    subjects=regexp(subjects,'listener[0-9]*.mat','match');
    subjects=subjects(~cellfun(@isempty,subjects));
    
    for si=1:length(subjects);
        others=subjects((1:length(subjects))~=si);
        for oi=1:length(others);
            load([expdir experiments{ei} '/fmri/mat/wholeBrain/' char(others{oi})]);
            
            if oi==1;
                data_others=zscore(data')';
            else;
                data_others=data_others+zscore(data')';
            end
        end
        clear data
        data=data_others/length(others);
        save([expdir experiments{ei} '/fmri/mat/wholeBrain/leave1out/zscore_leave_' char(subjects{si})  ],'data','-v7.3');
        clear data data_others
    end
    
end
