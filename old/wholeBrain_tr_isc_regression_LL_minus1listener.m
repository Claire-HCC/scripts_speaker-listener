
clear all;

tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

tic % 15 min

% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
voxind=111508;
lags=0;
for ei=3;%1:2;
    exp=experiments{ei};
    
    b=[];
    roi_ids=[];
    r2=[];
    F=[];
    
    clear data_mat

    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/' '/zscore_listenerAll.mat' ],'gdata','keptvox');
    
    %% isc or pattern
    % gdata=nanmean(gdata,1);
    vis=find(ismember(keptvox,voxind));
        gdata=gdata(vis,:,:);
        
    roi_voxn=1;
    roi_tn=size(gdata,2);
    keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags,0]));

    for vi=1:length(vis);
        for sp=1:size(gdata,3);
            listeners=1:size(gdata,3);
            listeners=listeners(listeners~=sp);
            
            for s=listeners;
                
                othersi=listeners(listeners~=s);
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
                
                [b(vi,:,sp,s),~,r(vi,keptT,sp,s),~,stats]=regress(y,X);
                
                r2(vi,:,sp,s)=stats(1);
                F(vi,:,sp,s)=stats(2);
                p(vi,:,sp,s)=stats(3);
                
                Y=X*b(vi,:,sp,s)';
                Y=reshape(Y,roi_voxn,length(keptT));
                y=reshape(y,roi_voxn,length(keptT));
                coupling(vi,:,sp,s)=zeros(1,roi_tn);
                coupling(vi,keptT,sp,s)=corr_col(y,Y);
                
                clear X
            end
            
        end
    end
    
    couplingz=0.5*log((1+coupling)./(1-coupling));
    save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/vox_regression_LLminus1L_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','r','F','r2','p','couplingz','lags','keptT','voxind','keptvox','gdata');
    clear b F p r2 rnames coupling
end

toc
