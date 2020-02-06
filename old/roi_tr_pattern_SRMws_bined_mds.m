clear all;
tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
rnames=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
rnames=strrep({rnames.name},'.mat','');

binSize=20;% tr;
binStep=1;
mds_dim=2;
type='corr';
col_speaker='r'; %'r'

for ei=3%1:2;
    exp=experiments{ei};
    
    ris=find(ismember(rnames,{'dPCC'}));
    % [~,ri]=max(herdm);
    
    for rii=1:length(ris);
        ri=ris(rii);
        rname=rnames{ri};
        fl=[expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '_listenersSRMws'  ];
        
        if exist([fl '.mat'])>0;
            load([expdir  exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '_listenersSRMws' ],'gdata');
            load([expdir  exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname ],'data');
            gdata=zscore(gdata,0,2);
            
            if size(gdata,1)>10;% & size(data2,1)>10;
                timeN=size(gdata2,2);
                subjN=size(gdata2,3);
                binN=length(1:binStep:(timeN-binSize+1));
                
                speaker_temp=[];
                gdata_temp=[];
                
                for bi=1:binN;
                    keptTi=((bi-1)*binStep+1):((bi-1)*binStep+binSize);
                    
                    temp=data2(:,keptTi);
                    speaker_temp(end+1,:)=temp(:);
                    for si=1:size(gdata2,3);
                        temp=gdata2(:,keptTi,si);
                        gdata_temp(end+1,:)=temp(:);
                    end
                end
                % temp=[speaker_temp ; gdata_temp];
                
                d=pdist(temp,type);
                Y=mdscale(double(d),mds_dim);
                
                % speaker_cord=Y(1:binN,:);
                listener_cord=Y;%((binN+1):end,:);
                listener_cord=reshape(listener_cord',mds_dim,subjN,binN);
                listener_cord=permute(listener_cord,[3 1 2]);
                
                save([expdir '/' exp '/fmri/' rname '_tr_pattern_SRMws_binSize' num2str(binSize) '_metric_mds_' type '_' col_speaker '_2d.mat'],'listener_cord','binN','binSize','binStep','rname')
                clear data2 data2 labels labels3
            end
        end
        
    end
end

