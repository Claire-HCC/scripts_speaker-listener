clear all;

tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
binSize=30; % tr;
binStep=1;

lags=0;
for ei=1:4;
    exp=experiments{ei};
    rnames=dir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/'  froidir '/zscore_listenerAll_*.mat']);
    rnames={rnames.name};
    rnames=strrep(rnames,'zscore_listenerAll_','');
    rnames=strrep(rnames,'.mat','');
    rnames=rnames';
    rnames=rnames(ismember(rnames,roi_table.region));
    
    b=[];
    roi_ids=[];
    r2=[];
    F=[];
    
    for ri=1:length(rnames);
        clear data_mat
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        
        keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags, 0])-binSize+1);
        binN=length(1:binStep:(length(keptT)-binSize+1));
                
        for sf=1:size(gdata,3);
            
             data=gdata(:,:,sf);
               listeners=1:size(gdata,3);
                listeners=listeners(listeners~=sf);
                g=nanmean(gdata(:,:,listeners),3);
                
                
            for bi=1:binN;
                keptTi=((bi-1)*binStep+1):((bi-1)*binStep+binSize);
                
                y=g;
                y=y(:,keptT(keptTi));
                y=y(:);
                
                for li=1:length(lags);
                    X(:,:,li)=data(:,keptT(keptTi)+lags(li));
                end
                
                X=reshape(X,roi_voxn*length(keptTi),length(lags));
                
                % centralized X
                X=X-mean(X);
                
                % add an constant
                X=[ones(size(X,1),1) X];
                
                [b(ri,:,sf),~,~,~,stats]=regress(y,X);
                
                r2(ri,keptT(keptTi(1)),sf)=stats(1);
                F(ri,keptT(keptTi(1)),sf)=stats(2);
                p(ri,keptT(keptTi(1)),sf)=stats(3);
                
                Y=X*b(ri,:,sf)';
                Y=reshape(Y,roi_voxn,length(keptTi));
                y=reshape(y,roi_voxn,length(keptTi));

                clear X
            end
        end
    end

    save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_SLfake_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','keptT','y','Y','rnames');
 %   clear b F p r2 rnames coupling
    
end
toc
