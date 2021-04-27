%б@всврв█бH
clear all;
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
rnames=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
rnames=strrep({rnames.name},'.mat','');

binSize=20;% tr;
binStep=1;

% type='corr';
% col_speaker='r'; %'r'

for ei=3%1:2;
    exp=experiments{ei};
    
    ris=find(ismember(rnames,{'dPCC'}));
    % [~,ri]=max(herdm);
    
    for rii=1:length(ris);
        ri=ris(rii);
        rname=rnames{ri};
        fl=[expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname  ];
        
        if exist([fl '.mat'])>0;
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname ],'data');
            
            if size(gdata,1)>10 & size(data,1)>10;
                timeN=size(gdata,2);
                subjN=size(gdata,3);
                binN=length(1:binStep:(timeN-binSize+1));
                
                speaker_temp=[];
                gdata_temp=[];
                
                for bi=1:binN;
                    keptTi=((bi-1)*binStep+1):((bi-1)*binStep+binSize);
                    
                    temp=data(:,keptTi);
                    speaker_temp(end+1,:)=temp(:);
                    for si=1:size(gdata,3);
                        temp=gdata(:,keptTi,si);
                        gdata_temp(end+1,:)=temp(:);
                    end
                end
                temp=[speaker_temp ; gdata_temp];
                
                  [coeff,score,latent,tsquared,explained,mu] = pca(temp);
                    speaker_score=score(1:binN,:);
                listener_score=score((binN+1):end,:);
                listener_score=reshape(listener_score',size(speaker_score,2),subjN,binN);
                listener_score=permute(listener_score,[3 1 2]);
            end
        end
    end
end
