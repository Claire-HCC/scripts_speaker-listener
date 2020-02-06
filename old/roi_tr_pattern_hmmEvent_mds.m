clear all;
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
rnames=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
rnames=strrep({rnames.name},'.mat','');
type='euclidean';
col_speaker='r'; %'r'
mds_dim=2;

K=28;
for ei=3%1:2;
    exp=experiments{ei};
    
    ris=find(ismember(rnames,{'dPCC'}));
    % [~,ri]=max(herdm);
    
    for rii=1:length(ris);
        ri=ris(rii);
        rname=rnames{ri};
        
        load([expdir '/scripts_speaker-listener/' exp '_' rname '_hmm_findSpeakerEventInListeners.mat'],'segments_S','segments_L');
        segments_S=segments_S(K,:);
        segments_L=segments_L(K,:,:);
        
        fl=[expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname   ];
        
        if exist([fl '.mat'])>0;
            load([expdir  exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname  ],'gdata');
            load([expdir  exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname ],'data');
            gdata=zscore(gdata,0,2);
            
            if size(gdata,1)>10 &  size(data,1)>10;
                timeN=size(gdata,2);
                subjN=size(gdata,3);
                evN=max(segments_S);
                
                speaker_temp=[];
                gdata_temp=[];
                
                for evi=1:evN;
                    speaker_temp(evi,:)=mean(data(:,segments_S==evi),2)';
                end
                
                for si=1:subjN;
                    for evi=1:evN;
                        gdata_temp(end+1,:)=squeeze(mean(gdata(:,segments_L(:,:,si)==evi,si),2));
                    end
                end
                temp=[speaker_temp ; gdata_temp];
                i=find(sum(isnan(temp),2)==0);
                temp2=temp(i,:);
                
                d=pdist(temp2,type);
                Y=mdscale(double(d),mds_dim);
                Y2=nan(size(temp,1),mds_dim);
                Y2(i,:)=Y;
                
                speaker_cord=Y2(1:evN,:);
                listener_cord=Y2((evN+1):end,:);
                listener_cord=reshape(listener_cord',mds_dim,evN,subjN);
                listener_cord=permute(listener_cord,[2 1 3]);
                
                save([expdir '/' exp '/fmri/' rname '_tr_pattern_hmmEvent_metric_mds_' type '_' col_speaker '_2d.mat'],'speaker_cord','listener_cord','evN','rname')
                clear data2 data2 labels labels3
            end
        end
        
    end
end

