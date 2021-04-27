clear all;
tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
rnames=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
rnames=strrep({rnames.name},'.mat','');

binSize=10;% tr;
binStep=1;

type='corr';
col_speaker='r'; %'r'
   mds_dim=1;
for ei=3:4%1:2;
    exp=experiments{ei};
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        
        if ismember(rname,{'pANG_L'}); %sherlock DLPFC_L
            
            fl=[expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname  ];
            
            if exist([fl '.mat'])>0;
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname ],'gdata');
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname ],'data');
                
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
                    
                    d=pdist(temp,type);
                 
                    % When starting with cmdscale, pieman vPCUN yielded
                    % this error message:
                    % Error using mdscale (line 288)
                    % Cannot initialize with 'cmdscale' when there are zero weights or missing data.
                    % Initialize with 'random' or provide explicit starting points.
                    
                    Y=mdscale(double(d),mds_dim,'start','random');
                    
                    speaker_cord=Y(1:binN,:);
                    listener_cord=Y((binN+1):end,:);
                    listener_cord=reshape(listener_cord',mds_dim,subjN,binN);
                    listener_cord=permute(listener_cord,[3 1 2]);
                    
                    save([expdir '/' exp '/fmri/animation/' rname '_pattern_binSize' num2str(binSize) '_mds_' type  '_' num2str(mds_dim)  'd.mat'],'speaker_cord','listener_cord','binN','binSize','binStep','rname','d')
                    
                end
            end
        end
    end
end

