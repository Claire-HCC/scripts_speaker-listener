clear all;

tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

tic % 15 min

% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;

lags=-6:-3;
voxind=111508;
for ei=3;%1:4;
    exp=experiments{ei};
    
    
    b=[];
    roi_ids=[];
    r2=[];
    F=[];
    
    clear data_mat
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_listenerAll.mat' ],'gdata','keptvox');
    vis=find(ismember(keptvox,voxind));
    gdata=gdata(vis,:,:);
    
    roi_voxn=1;
    roi_tn=size(gdata,2);
    keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags,0]));
    
    for vi=1:size(vis);
        for sf=1:size(gdata,3);
            data=gdata(vi,:,sf);
            
            listeners=1:size(gdata(vi,:,:),3);
            listeners=listeners(listeners~=sf);
            
            g=nanmean(gdata(vi,:,listeners),3);
            y=g;
            y=y(:,keptT);
            y=y(:);
            
            for li=1:length(lags);
                X(:,:,li)=data(:,keptT+lags(li));
            end
            
            X=reshape(X,roi_voxn*length(keptT),length(lags));
            
            % centralized X
            X=X-mean(X);
            
            % add an constant
            X=[ones(size(X,1),1) X];
            
            [b(vi,:,sf),~,r(vi,keptT,sf),~,stats]=regress(y,X);
            
            r2(vi,:,sf)=stats(1);
            F(vi,:,sf)=stats(2);
            p(vi,:,sf)=stats(3);
            
            Y=X*b(1,:,sf)';
            Y=reshape(Y,roi_voxn,length(keptT));
            y=reshape(y,roi_voxn,length(keptT));
            coupling(vi,:,sf)=zeros(1,roi_tn);
            coupling(vi,keptT,sf)=corr_col(y,Y);
            
            clear X
        end
    end
        
        couplingz=0.5*log((1+coupling)./(1-coupling));
        save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/vox_regression_SLfake_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r','r2','p','couplingz','lags','keptT','voxind','keptvox','y','Y');
        clear b F p r2 rnames coupling
        
    end
    toc
