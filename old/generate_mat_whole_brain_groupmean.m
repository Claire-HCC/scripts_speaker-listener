
% clear all
close all

loc='mypc';
set_parameters;


% %% generate the other images for self-other analysis
% for ei=1%%:2;
%
%     exp= experiments{ei};
%     % load([expdir exp '/fmri/mat/wholeBrain/listenerAll.mat']);
%     data_all=data;
%     clear data;
%
%     mask= (sum(sum(data_all==0,3),2))==0;
%     keptvox=find(mask==1);
%     data_all=data_all(keptvox,:,:);
%
%     subjn=size(data_all,3);
%     for subji=2:subjn;
%         others=find(~ismember(1:subjn,subji));
%         % do zscore twice will get the same results as zscore twice. So I
%         % jsut do it to be more conservative.
%         datam=mean(zscore(data_all(:,:,others),0,2),3);
%
%         data=zeros(voxn,size(datam,2));
%         data(keptvox,:)=datam;
%
%         fname=sprintf('%s/%s/fmri/mat/wholeBrain/leave1out/listener%02d_others.mat',expdir,exp,subji);
%         save(fname,'data');
%         clear data datam;
%     end
%     clear data_all data datam
% end


% for ei=1%%:2;

exp= experiments{ei};

subjects=dir([expdir experiments{ei} '/fmri/mat/wholeBrain/listener*mat']);
subjects={subjects.name}';
subjects=subjects(~cellfun(@isempty,(regexp(subjects,'[0-9]*'))));

subjn=length(subjects);

for si=1:subjn;
    
    subj=subjects{si};
    
    load([expdir exp '/fmri/mat/wholeBrain/' subj ]);
    data_all(:,:,si)=data;
end

data=data_all;
save([expdir experiments{ei} '/fmri/mat/wholeBrain/listenerAll.mat'],'data','-v7.3');

clear data

data_m=zeros(size(data_all(:,:,1)));
for si=1:size(data_all,3);
    % data_m=data_m+zscore(data_all(:,:,si)')';
    data_m=data_m+data_all(:,:,si);
end

data_m=data_m/size(data_all,3);
data=data_m;
save([expdir experiments{ei} '/fmri/mat/wholeBrain/listenerMean.mat'],'data');

clear data_all



%% divide subjects into two groups randomly and save the group averages;
% for ei=1%%:2;
%
%     exp= experiments{ei};
%
%     subjects=dir([expdir experiments{ei} '/fmri/mat/wholeBrain/listener*mat']);
%     subjects={subjects.name}';
%     subjects=subjects(~cellfun(@isempty,(regexp(subjects,'[0-9]*'))));
%
%     subjn=length(subjects);
%
%     [~,subj_permuted] = sort(rand(subjn,1));
%     subj_permuted=reshape(subj_permuted,2,subjn/2);
%
%     for gi=1:2;
%
%         for si=1:subjn/2;
%
%             subj=subjects{subj_permuted(gi,si)};
%
%             load([expdir exp '/fmri/mat/wholeBrain/' subj ]);
%             data_all(:,:,si)=data;
%
%             clear data
%         end
%
%         data_m=zeros(size(data_all(:,:,1)));
%
%         for si=1:size(data_all,3);
%             data_m=data_m+zscore(squeeze(data_all(:,:,si))')';
%         end
%
%         data_m=data_m/size(data_all,3);
%         data=data_m;
%         save([expdir experiments{ei} '/fmri/mat/wholeBrain/listenerZscoreMean_g' num2str(gi) '.mat'],'data');
%
%         clear data_all
%     end
% end