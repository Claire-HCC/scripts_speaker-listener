
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
timeBinLen=3;

lags=-4:4;
for ei=2:4;
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
        
        %% isc or pattern
        % gdata=nanmean(gdata,1);
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags,0]));
        
        for s=1:size(gdata,3);
            
            othersi=1:size(gdata,3);
            others=othersi(othersi~=s);
            others=nanmean(gdata(:,:,othersi),3);
            
            self=gdata(:,:,s);
            
            y=others;
            y=y(:,keptT);
            y=y(:);
            
            for li=1:length(lags);
                X(:,:,li)=self(:,keptT+lags(li));
            end
            
            X=reshape(X,roi_voxn*length(keptT),length(lags));
            
            % centralized X
            X=X-mean(X);
            
            % add an constant
            X=[ones(size(X,1),1) X];
            
            [b(ri,:,s),~,~,~,stats]=regress(y,X);
            
            r2(ri,:,s)=stats(1);
            F(ri,:,s)=stats(2);
            p(ri,:,s)=stats(3);
            
            Y=X*b(ri,:,s)';
            Y=reshape(Y,roi_voxn,length(keptT));
            y=reshape(y,roi_voxn,length(keptT));
            
            timeBinN=floor(size(Y,2)/timeBinLen);
            Y_bined=Y(:,1:timeBinN*timeBinLen);
            Y_bined=reshape(Y_bined(:),size(Y,1),timeBinLen,timeBinN);
            Y_bined=squeeze(mean(Y_bined,2));
            
            y_bined=y(:,1:timeBinN*timeBinLen);
            y_bined=reshape(y_bined(:),size(Y,1),timeBinLen,timeBinN);
            y_bined=squeeze(mean(y_bined,2));
            
            coupling(ri,:)=zeros(1,timeBinN);
            coupling(ri,1:timeBinN)=corr_col(y_bined,Y_bined);
            
            clear X
        end
        
    end
    couplingz=0.5*log((1+coupling)./(1-coupling));
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_bined'],'b','F','r2','p','couplingz','lags','rnames','keptT','timeBinN','timeBinLen');
    clear b F p r2 rnames coupling
end

toc
