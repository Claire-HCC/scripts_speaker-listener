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
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/' '/zscore_listenerAll.mat' ],'gdata');
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/' '/zscore_speaker.mat'],'data','keptvox');
    
    %% isc or pattern
    % gdata=nanmean(gdata,1);
    % data=nanmean(data,1);
    
    vis=find(ismember(keptvox,voxind));
        gdata=gdata(vis,:,:);
    data=data(vis,:);
    
    roi_voxn=1;
    roi_tn=size(gdata,2);

    
    for vi=1:size(data,1);
        for sp=1:size(gdata,3);
            listeners=1:size(gdata,3);
            listeners=listeners(listeners~=sp);
            g=nanmean(gdata(vi,:,listeners),3);
            
            y=g;
            keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags, 0]));
            y=y(:,keptT);
            y=y(:);
            
            for li=1:length(lags);
                X(:,:,li)=data(vi,keptT+lags(li));
            end
            
            X=reshape(X,roi_voxn*length(keptT),length(lags));
            
            % centralized X
            X=X-mean(X);
            
            % add an constant
            X=[ones(size(X,1),1) X];
            
            [b(vi,:,1,sp),~,r(vi,keptT,1,sp),~,stats]=regress(y,X);
            
            r2(vi,:,sp)=stats(1);
            F(vi,:,sp)=stats(2);
            p(vi,:,sp)=stats(3);
            
            Y=X*b(vi,:,sp)';
            Y=reshape(Y,roi_voxn,length(keptT));
            y=reshape(y,roi_voxn,length(keptT));
            coupling(vi,:,sp)=zeros(1,roi_tn);
            coupling(vi,keptT,sp)=corr_col(y,Y);
            
            clear X
        end
    end
    
    couplingz=0.5*log((1+coupling)./(1-coupling));
    save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/vox_regression_SLminus1L_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','r','F','r2','p','couplingz','lags','keptT','voxind','keptvox','Y','y');
    clear b F p r2 rnames coupling
end

toc
beep
