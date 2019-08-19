
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
[~,rnames_temp,~] = cellfun(@fileparts,cellstr(ls([expdir '/roi_mask/'  froidir '/mat/*.mat'])),'UniformOutput',0);
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

lags=-5:5;
for ei=1:2;
    exp=experiments{ei};
    
    b_SL=[];
    roi_ids=[];
    rii=1;
    for ri=1:length(rnames_temp);
        clear data_mat
        rname=rnames_temp{ri};
        fl=[expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname  ];
        
        if exist([fl '.mat'])>0;
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname ],'data');
            
            for s=1:size(gdata,3);
                
                self=gdata(:,:,s);
                othersi=1:size(gdata,3);
                othersi(othersi==s)=[];
                others=nanmean(gdata(:,:,othersi),3);
                
                keptT=(max(lags)+1):(size(data,2)-max(lags));
                others=others(:,keptT);
                
                for li=1:length(lags);
                    data_mat(:,:,li)=data(:,keptT+lags(li));
                end
                others=others(:);
                data_mat=reshape(data_mat,size(data_mat,1)*size(data_mat,2),length(lags));
                b_SL(rii,:,s)=ridge(others,data_mat,0);
                clear data_mat
            end
              rnames{rii,1}=rname;
            rii=rii+1;
          
        end
        
    end
    
    save([expdir '/' exp '/fmri/pattern_ridge/' timeUnit '/roi/' froidir '/beta_SL_' rname ],'b_SL','lags','rnames');
    % clear b_SL roi_ids rnames
end

toc
